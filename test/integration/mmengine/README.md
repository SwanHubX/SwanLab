# SwanLab integrate MMEngine test docs

## User`s guidance

There are two scripts to test intgration module

* `mmengine_train.py` for full test. It is a implement for cifar10 classification mission. Use resnet50 as our classifier and use mmengine with our training farmwork. Ref [mmengine docs](https://mmengine.readthedocs.io/en/latest/get_started/15_minutes.html).

* `mmengine_visualizer_import.py` is a simple nad efficient test scripts. It use `mmengine.registry` initialize mmengine 'visualizer' with swanlab backend.  Ref [mmegine docs](https://mmengine.readthedocs.io/en/latest/advanced_tutorials/visualization.html#customize-storage-backends-and-visualizers)

## Install

Refer to [mmengine official document](https://mmengine.readthedocs.io/en/latest/get_started/installation.html).

Install the environment with following command.

```sh
# with cuda12.1 or you can find torch version you want at pytorch.org
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121

pip install -U openmim
mim install mmengine
pip install swanlab
```

Or use `pip install -r requirements.txt` to install env. (hope it can work)

## Start Test

Simple test with init visualizer

```sh
python mmengine_visualizer_import.py
```

Full test with training mission

```sh
python mmengine_train.py
```