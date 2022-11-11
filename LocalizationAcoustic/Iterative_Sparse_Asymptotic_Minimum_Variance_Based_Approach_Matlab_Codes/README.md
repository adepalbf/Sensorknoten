# MATLAB Implementaton Code for SAMV 
Iterative Sparse Asymptotic Minimum Variance Based Approaches for Array Processing
Abeida, Habti, Qilin Zhang, Jian Li, and Nadjim Merabtine. “Iterative sparse asymptotic minimum variance based approaches for array processing.” IEEE Transactions on Signal Processing 61, no. 4 (2013): 933-944.

Paper official Link: http://dx.doi.org/10.1109/TSP.2012.2231676
Code Website: https://qilin-zhang.github.io/publications/
Paper PDF: https://qilin-zhang.github.io/_pages/pdfs/SAMVpaper.pdf?raw=true


If you find our code useful for your research, please cite
```
@article{abeida2013iterative,
  title={Iterative sparse asymptotic minimum variance based approaches for array processing},
  author={Abeida, Habti and Zhang, Qilin and Li, Jian and Merabtine, Nadjim},
  journal={Signal Processing, IEEE Transactions on},
  volume={61},
  number={4},
  pages={933--944},
  year={2013},
  publisher={IEEE}
}
```


### Where to start
plotDOA_est_MCtimes.m: it generates DOA estimation results in the paper. 

driver_parfor_MCs.m: It calls the function "Parfor_MC_SAMS.m".  "Parfor_MC_SAMS.m" takes advantage of an early version of MATLAB parallel computing toolbox (which changes the syntax from "matlabpool" to "parpool"), so some modification will be required to run "Parfor_MC_SAMS.m" in parallel with 8 cores. An alternative way to run it with a single core without parallel support is also feasible, just remove all "parallel related things" and replace "parfor" with the regular "for". 

### Contact
Please find the contact information at https://qilin-zhang.github.io/contact/
