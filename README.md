# fu_evapotranspiration
A python implementation of the model from:

Graamans, L., van den Dobbelsteen, A., Meinen, E., & Stanghellini, C. (2017). Plant factories; crop transpiration and energy balance. Agricultural Systems, 153, 138-147. [https://doi.org/10.1016/j.agsy.2017.01.003](https://doi.org/10.1016/j.agsy.2017.01.003)

The repository contains two mains files:

**evt_model.py**

This contains the model and can be executed in a terminal by cd'ing into the repository directory and running the command:

`python ./evt_model.py`

**test_evt_model.py**

This contains functions to test the code  in evt_model.py. The tests either test the functions against values given in the paper, or against the reference FAO Penman-Monteith crop transpiration model (on which the model si based): http://www.fao.org/3/X0490E/x0490e00.htm

The model is 95% complete, but currently cannot replicate the values from the paper - presumably due to confusion about the units being used.
