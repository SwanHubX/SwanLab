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
import logging
import swanlab
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence
from swanlab.intergration.utils.autologging import Response


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
        start_time: float,  # pass to comply with the protocol, but use response["created"] instead
        time_elapsed: float,
    ) -> Optional[Dict[str, Any]]:
        request = kwargs

        try:
            if response.get("object") == "edit":
                return self._resolve_edit(request, response, time_elapsed)
            elif response.get("object") == "text_completion":
                return self._resolve_completion(request, response, time_elapsed)
            elif response.get("object") == "chat.completion":
                return self._resolve_chat_completion(request, response, time_elapsed)
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
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Edit`."""
        request_str = f"\n\n**Instruction**: {request['instruction']}\n\n" f"**Input**: {request['input']}\n"
        choices = [f"\n\n**Edited**: {choice['text']}\n" for choice in response["choices"]]

        return self._resolve_metrics(
            request=request,
            response=response,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_completion(
        self,
        request: Dict[str, Any],
        response: Response,
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        request_str = f"\n\n**Prompt**: {request['prompt']}\n"
        choices = [f"\n\n**Completion**: {choice['text']}\n" for choice in response["choices"]]

        return self._resolve_metrics(
            request=request,
            response=response,
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_chat_completion(
        self,
        request: Dict[str, Any],
        response: Response,
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
            request_str=request_str,
            choices=choices,
            time_elapsed=time_elapsed,
        )

    def _resolve_metrics(
        self,
        request: Dict[str, Any],
        response: Response,
        request_str: str,
        choices: List[str],
        time_elapsed: float,
    ) -> Dict[str, Any]:
        """Resolves the request and response objects for `openai.Completion`."""
        results = [{"inputs": {"request": request_str}, "outputs": {"response": choice}} for choice in choices]
        metrics = self._get_metrics_to_log(request, response, results, time_elapsed)
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
        results: List[Any],
        time_elapsed: float,
    ) -> Metrics:
        model = response.get("model") or request.get("model")
        usage_metrics = self._get_usage_metrics(response, time_elapsed)

        usage = []
        for result in results:

            row = {
                "request": result.inputs["request"],
                "response": result.outputs["response"],
                "model": model,
                "start_time": datetime.datetime.fromtimestamp(response["created"]),
                "end_time": datetime.datetime.fromtimestamp(response["created"] + time_elapsed),
                "request_id": response.get("id", None),
                "api_type": response.get("api_type", "openai"),
            }
            result_text = swanlab.Text(row)

        metrics = Metrics(stats=result_text, usage=usage_metrics)
        return metrics

    @staticmethod
    def _convert_metrics_to_dict(metrics: Metrics) -> Dict[str, Any]:
        """Converts metrics to a dict."""
        metrics_dict = {
            "stats": metrics.stats,
            "trace": metrics.trace,
        }
        usage_stats = {f"usage/{k}": v for k, v in asdict(metrics.usage).items()}
        metrics_dict.update(usage_stats)
        return metrics_dict