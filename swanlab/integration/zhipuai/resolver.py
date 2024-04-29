#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-26 17:07:44
@File: swanlab\integration\zhipuai\resolver.py
@IDE: vscode
@Description:
    ZhipuAI Resolver
"""
import datetime
import io
import os
import sys
import swanlab
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence

current_dir = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.join(current_dir, "..")
sys.path.append(relative_path)

from integration_utils.autologging import Response


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


class ZhipuAIClientResponseResolver:
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

            response = self.chat_completion_to_dict(response)
            return self._resolve_chat_completion(request, response, lib_version, time_elapsed)

        except Exception as e:
            print(f"Failed to resolve request/response: {e}")
        return None

    def _resolve_chat_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        lib_version: str,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `zhipuai.ZhipuAI().Completion`."""
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
        """Resolves the request and response objects for `zhipuai.ZhipuAI().Completion`."""
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
        model = response.get("model")
        usage_metrics = self._get_usage_metrics(response, time_elapsed)

        usage = []

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

    def chat_completion_to_dict(self, completion):
        data = {
            "model": completion.model,
            "created": completion.created,
            "choices": [
                {
                    "index": choice.index,
                    "finish_reason": choice.finish_reason,
                    "message": {
                        "content": choice.message.content,
                        "role": choice.message.role,
                        "tool_calls": choice.message.tool_calls,
                    },
                }
                for choice in completion.choices
            ],
            "request_id": completion.request_id,
            "id": completion.id,
            "usage": {
                "prompt_tokens": completion.usage.prompt_tokens,
                "completion_tokens": completion.usage.completion_tokens,
                "total_tokens": completion.usage.total_tokens,
            },
        }
        return data
