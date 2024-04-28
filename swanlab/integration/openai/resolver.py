#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-23 11:25:22
@File: swanlab\intergration\openai\resolver.py
@IDE: vscode
@Description:
    OpenAI Resolver
"""
import datetime
import io
import swanlab
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence
from swanlab.integration.integration_utils.autologging import Response


@dataclass
class UsageMetrics:
    elapsed_time: float = None
    prompt_tokens: int = None
    completion_tokens: int = None
    total_tokens: int = None


@dataclass
class Metrics:
    usage: UsageMetrics = None
    stats: swanlab.Text = None


class OpenAIRequestResponseResolver:
    def __init__(self):
        self.define_metrics_called = False

    def __call__(
        self,
        args: Sequence[Any],
        kwargs: Dict[str, Any],
        response: Response,
        lib_version: str,
        start_time: float,  # pass to comply with the protocol, but use response["created"] instead
        time_elapsed: float,
    ) -> Optional[Dict[str, Any]]:
        request = kwargs
        lib_version = lib_version
        try:
            if response.get("object") == "edit":
                return self._resolve_edit(request, response, lib_version, time_elapsed)
            elif response.get("object") == "text_completion":
                return self._resolve_completion(request, response, lib_version, time_elapsed)
            elif response.get("object") == "chat.completion":
                return self._resolve_chat_completion(request, response, lib_version, time_elapsed)
            else:
                # todo: properly treat failed requests
                print(f"Unsupported OpenAI response object: {response.get('object')}")
        except Exception as e:
            print(f"Failed to resolve request/response: {e}")
        return None

    def _resolve_edit(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Edit`."""
        request_str = f"\n\n**Instruction**: {request['instruction']}\n\n" f"**Input**: {request['input']}\n"
        choices = [f"\n\n**Edited**: {choice['text']}\n" for choice in response["choices"]]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        request_str = f"\n\n**Prompt**: {request['prompt']}\n"
        choices = [f"\n\n**Completion**: {choice['text']}\n" for choice in response["choices"]]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_chat_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        prompt = io.StringIO()
        for message in request["messages"]:
            prompt.write(f"\n\n**{message['role']}**: {message['content']}\n")
        request_str = prompt.getvalue()

        choices = [
            f"\n\n**{choice['message']['role']}**: {choice['message']['content']}\n" for choice in response["choices"]
        ]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_metrics(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        request_str: str,
        choices: List[str],
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        results = [{"inputs": {"request": request_str}, "outputs": {"response": choice}} for choice in choices]
        metrics = self._get_metrics_to_log(request, response, lib_version, results, time_elapsed)
        return self._convert_metrics_to_dict(metrics)

    @staticmethod
    def _get_usage_metrics(response: Response, time_elapsed: float) -> UsageMetrics:
        """Gets the usage stats from the response object."""
        if response.get("usage"):
            usage_stats = UsageMetrics(**response["usage"])
        else:
            usage_stats = UsageMetrics()
        usage_stats.elapsed_time = time_elapsed
        return usage_stats

    def _get_metrics_to_log(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        results: List[Any],
        time_elapsed: float,
    ) -> Metrics:
        model = response.get("model") or request.get("model")
        usage_metrics = self._get_usage_metrics(response, time_elapsed)

        usage = []

        for result in results:
            swanlab.log(
                {
                    "result": [
                        swanlab.Text(caption="request", data=str(request["messages"])),
                        swanlab.Text(caption="response", data=response.choices[0].message.content),
                        swanlab.Text(caption="lib_version", data=lib_version),
                        swanlab.Text(caption="request_id", data=response.id),
                        swanlab.Text(caption="model", data=response.model),
                        swanlab.Text(caption="completion_tokens", data=response.usage.completion_tokens),
                        swanlab.Text(caption="prompt_tokens", data=response.usage.prompt_tokens),
                        swanlab.Text(caption="total_tokens", data=response.usage.total_tokens),
                        swanlab.Text(
                            caption="start_time", data=str(datetime.datetime.fromtimestamp(response["created"]))
                        ),
                        swanlab.Text(
                            caption="end_time",
                            data=str(datetime.datetime.fromtimestamp(response["created"] + time_elapsed)),
                        ),
                        swanlab.Text(caption="used_time", data=time_elapsed),
                    ]
                }
            )

            row = {
                "request": str(request["messages"]),
                "response": response.choices[0].message.content,
                "lib_version": lib_version,
                "model": model,
                "start_time": datetime.datetime.fromtimestamp(response["created"]),
                "end_time": datetime.datetime.fromtimestamp(response["created"] + time_elapsed),
                "request_id": response.get("id", None),
                "api_type": response.get("api_type", "openai"),
            }
        result_text = swanlab.Text(caption="result_row_dict", data=str(row))

        metrics = Metrics(stats=result_text, usage=usage_metrics)
        return metrics

    @staticmethod
    def _convert_metrics_to_dict(metrics: Metrics) -> Dict[str, Any]:
        """Converts metrics to a dict."""
        metrics_dict = {
            "stats": metrics.stats,
        }
        usage_stats = {f"usage/{k}": v for k, v in asdict(metrics.usage).items()}
        metrics_dict.update(usage_stats)
        return metrics_dict


class OpenAIClientResponseResolver:
    def __init__(self):
        self.define_metrics_called = False

    def __call__(
        self,
        args: Sequence[Any],
        kwargs: Dict[str, Any],
        response: Response,
        lib_version: str,
        start_time: float,  # pass to comply with the protocol, but use response["created"] instead
        time_elapsed: float,
    ) -> Optional[Dict[str, Any]]:
        request = kwargs
        lib_version = lib_version

        try:
            if response.object == "edit":
                return self._resolve_edit(request, response, lib_version, time_elapsed)

            elif response.object == "text_completion":
                response = self.completion_to_dict(response)
                return self._resolve_completion(request, response, lib_version, time_elapsed)
            elif response.object == "chat.completion":
                response = self.chat_completion_to_dict(response)
                return self._resolve_chat_completion(request, response, lib_version, time_elapsed)
            else:
                # todo: properly treat failed requests
                print(f"Unsupported OpenAI response object: {response.get('object')}")
        except Exception as e:
            print(f"Failed to resolve request/response: {e}")
        return None

    def _resolve_edit(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Edit`."""
        request_str = f"\n\n**Instruction**: {request['instruction']}\n\n" f"**Input**: {request['input']}\n"
        choices = [f"\n\n**Edited**: {choice['text']}\n" for choice in response["choices"]]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.OpenAI().Completion`."""
        request_str = f"\n\n**Prompt**: {request['prompt']}\n"
        response_choices = response.get("choices")
        choices = [f"\n\n**Completion**: {choice['text']}\n" for choice in response_choices]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_chat_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.OpenAI().Completion`."""
        prompt = io.StringIO()
        for message in request["messages"]:
            prompt.write(f"\n\n**{message['role']}**: {message['content']}\n")
        request_str = prompt.getvalue()
        response_choices = response.get("choices")
        choices = [
            f"\n\n**{choice['message']['role']}**: {choice['message']['content']}\n" for choice in response_choices
        ]

        return self._resolve_metrics(
            request=request,
            response=response,
            lib_version=lib_version,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_metrics(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        request_str: str,
        choices: List[str],
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        results = [{"inputs": {"request": request_str}, "outputs": {"response": choice}} for choice in choices]
        metrics = self._get_metrics_to_log(request, response, lib_version, results, time_elapsed)
        return self._convert_metrics_to_dict(metrics)

    @staticmethod
    def _get_usage_metrics(response: Response, time_elapsed: float) -> UsageMetrics:
        """Gets the usage stats from the response object."""
        if response.get("usage"):
            usage_stats = UsageMetrics(**response["usage"])
        else:
            usage_stats = UsageMetrics()
        usage_stats.elapsed_time = time_elapsed
        return usage_stats

    def _get_metrics_to_log(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        results: List[Any],
        time_elapsed: float,
    ) -> Metrics:
        model = response.get("model") or request.get("model")
        usage_metrics = self._get_usage_metrics(response, time_elapsed)

        usage = []

        if response["object"] == "chat.completion":
            for result in results:
                swanlab.log(
                    {
                        "result": [
                            swanlab.Text(caption="request", data=str(request["messages"])),
                            swanlab.Text(caption="response", data=response["choices"][0]["message"]["content"]),
                            swanlab.Text(caption="lib_version", data=lib_version),
                            swanlab.Text(caption="request_id", data=response["id"]),
                            swanlab.Text(caption="model", data=response["model"]),
                            swanlab.Text(caption="completion_tokens", data=response["usage"]["completion_tokens"]),
                            swanlab.Text(caption="prompt_tokens", data=response["usage"]["prompt_tokens"]),
                            swanlab.Text(caption="total_tokens", data=response["usage"]["total_tokens"]),
                            swanlab.Text(
                                caption="start_time", data=str(datetime.datetime.fromtimestamp(response["created"]))
                            ),
                            swanlab.Text(
                                caption="end_time",
                                data=str(datetime.datetime.fromtimestamp(response["created"] + time_elapsed)),
                            ),
                            swanlab.Text(caption="used_time", data=time_elapsed),
                        ]
                    }
                )

                row = {
                    "request": str(request["messages"]),
                    "response": response["choices"][0]["message"]["content"],
                    "lib_version": lib_version,
                    "model": model,
                    "start_time": datetime.datetime.fromtimestamp(response["created"]),
                    "end_time": datetime.datetime.fromtimestamp(response["created"] + time_elapsed),
                    "request_id": response.get("id", None),
                    "api_type": response.get("api_type", "openai"),
                }

        # elif response["object"] == "text_completion":
        else:
            for result in results:
                swanlab.log(
                    {
                        "result": [
                            swanlab.Text(caption="request", data=str(request["prompt"])),
                            swanlab.Text(caption="response", data=response["choices"][0]["text"]),
                            swanlab.Text(caption="lib_version", data=lib_version),
                            swanlab.Text(caption="request_id", data=response["id"]),
                            swanlab.Text(caption="model", data=response["model"]),
                            swanlab.Text(caption="completion_tokens", data=response["usage"]["completion_tokens"]),
                            swanlab.Text(caption="prompt_tokens", data=response["usage"]["prompt_tokens"]),
                            swanlab.Text(caption="total_tokens", data=response["usage"]["total_tokens"]),
                            swanlab.Text(
                                caption="start_time", data=str(datetime.datetime.fromtimestamp(response["created"]))
                            ),
                            swanlab.Text(
                                caption="end_time",
                                data=str(datetime.datetime.fromtimestamp(response["created"] + time_elapsed)),
                            ),
                            swanlab.Text(caption="used_time", data=time_elapsed),
                        ]
                    }
                )

                row = {
                    "request": str(request["prompt"]),
                    "response": response["choices"][0]["text"],
                    "lib_version": lib_version,
                    "model": model,
                    "start_time": datetime.datetime.fromtimestamp(response["created"]),
                    "end_time": datetime.datetime.fromtimestamp(response["created"] + time_elapsed),
                    "request_id": response.get("id", None),
                    "api_type": response.get("api_type", "openai"),
                }

        result_text = swanlab.Text(caption="result_row_dict", data=str(row))

        metrics = Metrics(stats=result_text, usage=usage_metrics)
        return metrics

    @staticmethod
    def _convert_metrics_to_dict(metrics: Metrics) -> Dict[str, Any]:
        """Converts metrics to a dict."""
        metrics_dict = {
            "stats": metrics.stats,
        }
        usage_stats = {f"usage/{k}": v for k, v in asdict(metrics.usage).items()}
        metrics_dict.update(usage_stats)
        return metrics_dict

    def chat_completion_to_dict(self, chat_completion):
        chat_completion_dict = {
            "id": chat_completion.id,
            "choices": [
                {
                    "finish_reason": choice.finish_reason,
                    "index": choice.index,
                    "logprobs": choice.logprobs,
                    "message": {
                        "content": choice.message.content,
                        "role": choice.message.role,
                        "function_call": choice.message.function_call,
                        "tool_calls": choice.message.tool_calls,
                    },
                }
                for choice in chat_completion.choices
            ],
            "created": chat_completion.created,
            "model": chat_completion.model,
            "object": chat_completion.object,
            "system_fingerprint": chat_completion.system_fingerprint,
            "usage": {
                "completion_tokens": chat_completion.usage.completion_tokens,
                "prompt_tokens": chat_completion.usage.prompt_tokens,
                "total_tokens": chat_completion.usage.total_tokens,
            },
        }
        return chat_completion_dict

    def completion_to_dict(self, completion):
        completion_dict = {
            "id": completion.id,
            "choices": [
                {
                    "finish_reason": choice.finish_reason,
                    "index": choice.index,
                    "logprobs": choice.logprobs,
                    "text": choice.text,
                }
                for choice in completion.choices
            ],
            "created": completion.created,
            "model": completion.model,
            "object": completion.object,
            "system_fingerprint": completion.system_fingerprint,
            "usage": {
                "completion_tokens": completion.usage.completion_tokens,
                "prompt_tokens": completion.usage.prompt_tokens,
                "total_tokens": completion.usage.total_tokens,
            },
        }
        return completion_dict
