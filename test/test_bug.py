import swanlab

# Initialize
swanlab.init()

for epoch in range(1, 20):
    print("epoch", epoch, "|")
    # Tracking index: `epoch`
    swanlab.log({"epoch": epoch})

print("end")
print("===========", end="")
print("--------------------", end="")
