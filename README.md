# SolutionQuartic
Matlab code for quadratic equations in two variables. 求解二元二次方程组的Matlab代码

Much faster than Matlab vpasolve, useful to complexity and running time evaluation. 求解速度远快于Matlab vpasolve, 可用于评估计算复杂度和实际运行时间


Use the following code to call this function 使用以下代码调用该函数:

A1=[a1,b1,c1,d1,e1,f1]; % coefficient of the 1st quadratic equation 第一个二次方程的系数 a1·x^2 + b1·xy + c1·^2 + d1·x + e1·y + f1=0

B1=[a2,b2,c2,d2,e2,f2]; % coefficient of the 2nd quadratic equation 第二个二次方程的系数 a2·x^2 + b2·xy + c2·y^2 + d2·x + e2·y + f2=0

[x,y]=solvequartic(A1,B1);


**If you use or are inspired by this code, please give credit to the following paper. 如果使用本代码或受到启发，请引用以下论文**

**S. Zhao, X.-P. Zhang, X. Cui, and M. Lu, “A closed-form localization method utilizing pseudorange measurements from two non-synchronized positioning systems,” IEEE Internet Things J., vol. 8, no. 2, pp. 1082–1094, Jan. 2021.**

BibTex form:

@article{zhao2020closed,
  title={A Closed-form Localization Method Utilizing Pseudorange Measurements from Two Non-synchronized Positioning Systems},
  author={Zhao, Sihao and Zhang, Xiao-Ping and Cui, Xiaowei and Lu, Mingquan},
  journal={IEEE Internet of Things Journal},
  year={2020},
  publisher={IEEE}
}
