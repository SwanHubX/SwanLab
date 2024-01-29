from swanlab.db import Tag

# print(type(Tag.create_tag(1, "ikunskun-21")))

tags = Tag.get_tag(1)

print(tags)
