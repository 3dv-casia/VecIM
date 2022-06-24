# Vectorized Indoor Surface Reconstruction from 3D Point Cloud with Multistep 2D Optimization

## Introduction

This is the implementation of **VecIM** (ISPRS JPRS 2021), a multistep and versatile indoor LoD2 vectorized reconstruction pipeline without the Manhattan or Atlanta world assumptions. The core idea behind this method is the combination of a sequence of 2D segment or 
cell assembly problems defined as global optimizations while reducing the reconstruction complexity. For more details, please refer to the following [paper](https://www.sciencedirect.com/science/article/abs/pii/S0924271621001222):

    Jiali Han, Mengqi Rong, Hanqing Jiang, Hongmin Liu, Shuhan Shen,
    Vectorized indoor surface reconstruction from 3D point cloud with multistep 2D optimization,
    ISPRS Journal of Photogrammetry and Remote Sensing, 177: 57-74, 2021.

![image](https://github.com/ShuhanShen/VecIM/blob/main/images/pipeline.PNG)

## Run VecIM
### Dependencies
- CGAL (v4.11 has been tested)
- OPENCV (v3.4.3 has been tested)
- GLOG
- MPRF
- GMP

### Build
The project is built on CMake, and there is more than one way to compile it. The following is an example on Linux or macOS system:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_BUILD_TYPE=Release ..
    $ make

### Run
For a certain scene, the permanent structures including facade, floor, ceiling and cylinder (if any) are first segmented from the whole point cloud in the pre-processing. **VecIM** takes the point cloud of different components as the input and outputs a LoD2 indoor vectorized model.

Put different semantic point cloud (**facade.ply, floor.ply, ceiling.ply, cylinder.ply**) and parameter configuration file in the same folder with the path PATH, and then run:

    ./Facade/Facade PATH

An example of parameter configuration file is provided named as **config_modeling.xml**.

## License
This work is licensed under the [GNU General Public License v3.0](https://github.com/ShuhanShen/VecIM/blob/main/LICENSE).

## Acknowledgement
Part of the project is inspired by [PolyFit](https://github.com/LiangliangNan/PolyFit).

## Citation
If you use **VecIM** in a scientific work, please consider citing:

    @article{HAN2021vecim,
    title={Vectorized indoor surface reconstruction from 3D point cloud with multistep 2D optimization},
    author={Han, Jiali and Rong, Mengqi and Jiang, Hanqing and Liu, Hongmin and Shen, Shuhan},
    journal={ISPRS Journal of Photogrammetry and Remote Sensing},
    volume={177},
    pages={57--74},
    year={2021}
    }
