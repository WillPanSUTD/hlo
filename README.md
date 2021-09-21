# hlo
This is the source code for paper: Pan W, Lu X, Gong Y, et al. HLO: Half-kernel Laplacian operator for surface smoothing[J]. Computer-Aided Design, 2020, 121: 102807.

### Usage

Directly use the 'HLO.exe' file.


put the noisy mesh file under the same folder (e.g., bunny.off, supported format: obj, off, stl, wrl, ply, mesh)


run: 

```
HLO.exe bunny.off 5
```

'bunny.off' is the input mesh, 5 is the number of vertex update iteration.
You will see the denoised mesh name by 'xxx_denoised.obj' inside the same folder. 


If you find this code is useful for your research, please cite it by

@article{pan2020hlo,
  title={HLO: Half-kernel Laplacian operator for surface smoothing},
  author={Pan, Wei and Lu, Xuequan and Gong, Yuanhao and Tang, Wenming and Liu, Jun and He, Ying and Qiu, Guoping},
  journal={Computer-Aided Design},
  volume={121},
  pages={102807},
  year={2020},
  publisher={Elsevier}
}
