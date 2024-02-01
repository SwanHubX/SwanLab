from swanlab.db import Chart, Source, Tag, Display
from swanlab.db import connect


connect()

# Chart.create("hello", "world", experiment_id=6, type="audio", reference="time")
# Tag.create(6, "hello", "audio")
# Source.create(tag_id=13, chart_id=13)
Display.create(chart_id=13, namespace_id=6)
