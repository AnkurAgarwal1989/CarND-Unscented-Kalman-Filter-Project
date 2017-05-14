# Unscented Kalman Filter Project

UKF equations and code. Self-Driving Car Engineer Nanodegree Program

---

## Dependencies

* Eigen
* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`

## UKF Intuition:

We are tracking an object in 2D space. We are following the CTRV (Constant Turn Radius and Velocity) model as it seems to fit vehicular motion. Most of the time vehicles have a *constant velocity* and when turning also have a *constant turn radius*. To track such behavior, we need:

 - 2D position, $p_x$ and $p_y$.
 - a velocity, $v$.
 - a (yaw) heading angle, $\psi$.
 - rate of change of heading (yaw_rate), $\dot{\psi}$

Assembling into a state vector, $\mathbf{x} =  [p_x \,\, p_y \,\, v \,\, \psi \,\, \dot{\psi}]$



To "accomodate" changing velocity and changing yaw rate, we will add these as noise to the state.
$$
\text{acceleration}, \space\space a = \dot{v}\\
\text{yaw_rate_rate}, \space\space yaw_{dd} = \ddot{\psi}
$$

Unlike EKF which works by linearizing a non-linear system using a Taylor Approximation, UKF solves the problem by using special points that represent the function. These points are called ***Sigma Points***. We transform these points by the non-linear function and then recreate a Gaussian distribution.



### State Math

- Known: CTRV model


- State, $\underset{n_{x}\times 1} {\mathbf{x}} =  [p_x \,\, p_y \,\, v \,\, \psi \,\, \dot{\psi}]$

- We wish to find the function $f$, such that: $\mathbf{x_{k+1}}=f(\mathbf{x_k, \nu_k})$, where $\nu$ is the noise and $\mathbf{_k}$ is the time step. We can do so by calculating how the function changes w.r.t. time and then make predictions with $\Delta{t}$. Calculating $\dot{\mathbf{x}}$:

  - $\dot{p_x} = v_x = \cos{\psi}.v$ , Constant velocity and constant turning radius.
  - $\dot{p_y} = v_y = \sin{\psi}.v$ , Constant velocity and constant turning radius.
  - $\dot{v} = 0$ , Constant velocity.
  - $\dot{\psi} = \dot{\psi}$ , Already a part of the state.
  - $\ddot{\psi} = 0$ , Constant turning radius.

- Now, $\mathbf{x_{k+1}}=\mathbf{x_{k+1}} + \int_{t_k}^{t_{k+1}} \dot{\mathbf{x}}(t) dt$

  - **Calculating integrals for State:**

    - $\int_{t_k}^{t_{k+1}}\dot{\mathbf{p_x}}(t)dt$

    - $$
      \int_{t_k}^{t_{k+1}}\dot{\mathbf{p_x}}(t)dt = \int_{t_k}^{t_{k+1}}\cos{\psi_k}(t).v_k(t)dt\\
      = v_k\int_{t_k}^{t_{k+1}}\cos{\psi_k}(t)dt  \,\,\,\, ,v_k\, \mbox{remains constant so can be taken out}\\
      = v_k\int_{t_k}^{t_{k+1}}\cos({\psi_k + \dot{\psi_k}.(t-t_k)})dt  \,\,\,\, , \mbox{$\dot{\psi}$ is a known constant}\\
      \mbox{Substituting ${\psi_k + \dot{\psi_k}.(t-t_k)}$ with $x$}, \,\,\,\,\,\, dt = \frac{dx}{\dot{\psi_{k}}}, \,\,\,\,\,\ \mbox{lower limit}=\psi_k, \,\,\,\,\,\, \mbox{upper limit}=\psi_k+\dot{\psi_k}\Delta{t}\\
      = \frac{v_k}{\dot{\psi_k}}\int_{\psi_k}^{\psi_k+\dot{\psi_k}\Delta{t}} \cos(x)dx  \\
      =\fbox{$ \frac{v_k}{\dot{\psi_k}} (\sin(\psi_k+\dot{\psi_k}\Delta{t}) - \sin(\psi_k)  )$}\\
      \text{Additionally, to deal with $\dot{\psi_k} = 0$, we need to model system behavior for straight line motion}\\
      \fbox{$ \cos{\psi_k}.v_k\Delta{t}$} \,\,\,\,\,\,\,\, \text{when }\dot{\psi}=0\\
      $$

    - $\int_{t_k}^{t_{k+1}}\dot{\mathbf{p_y}}(t)dt$
      $$
      \text{Integration via substitution}\\
      \fbox{$ \frac{v_k}{\dot{\psi_k}} (-\cos(\psi_k+\dot{\psi_k}\Delta{t}) + \cos(\psi_k)  )$}\\
      \fbox{$ \sin{\psi_k}.v_k\Delta{t}$} \,\,\,\,\,\,\,\, \text{when }\dot{\psi}=0\\
      $$

    - $\int_{t_k}^{t_{k+1}}\dot{v}(t)dt$
      $$
      \text{Constant Velocity}\\
      \fbox{$\,\,\,\,\,$0$\,\,\,\,\,$}\\
      $$
      ​

    - $\int_{t_k}^{t_{k+1}}\dot{\psi}(t)dt$
      $$
      \text{Constant Turning Radius}\\
      =\dot{\psi_k}\int_{t_k}^{t_{k+1}}dt  \\
      \fbox{$\dot{\psi_k}\Delta{t}$}\\
      $$
      ​

    - $\int_{t_k}^{t_{k+1}}\ddot{\psi}(t)dt$
      $$
      \text{Constant Turning Radius}\\
      \fbox{$\,\,\,\,\,$0$\,\,\,\,\,$}\\
      $$

  - **Calculating integrals for Noise:**

    - Effect on $ \mathbf{p_x}(t): \int_{t_k}^{t_{k+1}}\nu_{\dot{v},k}dt$

    - Effect on $  \mathbf{p_y}(t): \int_{t_k}^{t_{k+1}}\nu_{\dot{v},k}dt$

    - Effect on $ v(t): \int_{t_k}^{t_{k+1}}\nu_{\dot{v},k}dt$
      $$
      \fbox{$\,\,\,\,\, \nu_{\dot{v},k}.\Delta{t} \,\,\,\,\,$}\\
      $$
      ​

    - Effect on $\psi(t) : \iint_{t_k}^{t_{k+1}}\nu_{\ddot{\psi}, k}dt$
      $$
      \fbox{$\frac{1}{2} \nu_{\ddot{\psi},k}.\Delta{t}^2$}\\
      $$
      ​

    - Effect on $ \dot{\psi}(t) : \int_{t_k}^{t_{k+1}}\nu_{\ddot{\psi}, k}dt$
      $$
      \fbox{$\,\,\,\,\, \nu_{\ddot{\psi},k}.\Delta{t} \,\,\,\,\,$}\\
      $$



### Steps

#### Predict

- Augment State and Process Matrix:
  - $\underset{n_{x}\times 1} {\mathbf{x}} =  [p_x \,\, p_y \,\, v \,\, \psi \,\, \dot{\psi}] \to \underset{n_{aug}\times 1} {\mathbf{x}} =  [p_x \,\, p_y \,\, v \,\, \psi \,\, \dot{\psi} \,\, \dot{v}  \,\,  \ddot{\psi} ]$
  - Process matrix, $\underset{n_x\times n_x} {\mathbf{P}}$ becomes $\underset{n_{aug}\times n_{aug}} {\begin{bmatrix} \mathbf{P} & 0\\ 0 & \mathbf{Q}\end{bmatrix}}$, where $\mathbf{Q}=\begin {bmatrix} \sigma^2_{} \dot{v}& 0 \\ 0 & \sigma^2_{ \ddot{\psi}}\end{bmatrix}$

- Generate a set of Sigma Points $\mathbf{\chi}$:

  - The spreading factor, $\lambda$ is a hyper-parameter that can be tuned. It controls how far from the mean our calculated sigma points will be. Larger the $\lambda$, larger the spread.

  - $\mathbf{\chi}$ is of size $n_{aug} \times (2n_{aug}+1)$. We form the points by referring to 

  - $$
    \mathbf{\chi} = \begin{array} {c|c} \mathbf{x} \space\space & \space\space \mathbf{x}+\sqrt{(\lambda+n_{aug})\mathbf{P}} \space\space&\space\space \mathbf{x}-\sqrt{(\lambda+n_{aug})\mathbf{P}} \end{array}
    $$

  - $\mathbf{\chi}^{[0]} = \mathbf{x}$

  - $\mathbf{\chi}^{[i+1]} = \mathbf{x}+\sqrt{(\lambda+n_{aug})\mathbf{P}^{[i]}} \,\,\,\,\,\,\, \mbox{ for }  i=0,...,n-1$

  - $\mathbf{\chi}^{[i+n_{aug}+1]} = \mathbf{x}+\sqrt{(\lambda+n_{aug})\mathbf{P}^{[i]}}  \,\,\,\,\,\,\, \mbox{ for }  i=0,...,n-1$

- Calculate weights for the Sigma Points:

  - Each of the $n_{aug} \times (2n_{aug}+1)$ sigma points has a weight attached to it, depending on it's distance from the mean of the distribution.
  - $\omega^{[0]}=\frac{\lambda}{\lambda+n_{aug}}$
  - $\omega^{[i]}=\frac{1}{2(\lambda+n_{aug})} \,\,\,\,\, \mbox{ for } i=1,...,2n$

- Make predictions by transforming the sigma points according to the state equations: