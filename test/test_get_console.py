from io import StringIO
import sys


class ConsoleCapture:
    def __enter__(self):
        self.original_stdout = sys.stdout
        sys.stdout = self._stdout = StringIO()
        return self

    def __exit__(self, *args):
        sys.stdout = self.original_stdout

    def get_captured_text(self):
        return self._stdout.getvalue()


# 用法示例
with ConsoleCapture() as cc:
    print("Hello, this will be captured.")
    print("So will this.")

captured_text = cc.get_captured_text()
print("Captured text:")
print(captured_text)
