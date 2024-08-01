# Modeling, Simulation, and Experimental Testing of a Single-legged Hopper
The following repository contains the Julia and Simulink files used in developing a physical single-legged forward hopper. The details in developing stable hopping for the single leg begins with modeling the dynamics, then simulation to develop a controller, and finally experimental testing to implement what was learned from the simulation and modeling to achieve stable forward hopping. This repository contains all the code developed for this project.

![Testing Sim2](Simulation/Animations/Bounding_Simulation2.gif)

## Modeling
A generic modeling method utilizing the Euler-Lagrange method derives the equations of motion in julia and saves the results to text files. This method is used to derive a variety of energy based systems. To establishe a symbolic mathematical model for the single leg's dynamics, The model of the leg is defined in two parts: the dynamics of flight and the dynamics of stance. 

Links to Julia derevation code.
- [Modeling Flight](Modeling/Modeling%20Software/SymPy%20DP%20Derivation2.0.jl)
- [Modeling Stance](Modeling/Modeling%20Software/SymPy%20DP%20Ground%20wSensor%20Derivation3.0.jl)
- [Modeling Inverted Pendulum](Modeling/Modeling%20Software/Inverted%20Pendulum%20Derivation.jl)
- [Modeling Triple Inverted Pendulum](Modeling/Modeling%20Software/SymPy%20Triple%20Pendulum%20Derivation.jl)

Links to derevation text files. The Equasions are written to be easily copied and pasted into code.
- [Flight Equasion](Modeling/Model%20Equasions/Flight.txt) for Julia
- [Stance Equasion](Modeling/Model%20Equasions/Stance.txt) for Julia
- Inverted Pendulum Equasion
  - [Generic Form](Modeling/Model%20Equasions/InvertedPendulum1.txt)
  - [For Matlab](Modeling/Model%20Equasions/InvertedPendulum2.txt)
- Triple Inverted Pendulum Equasion
  - [Generic Form](Modeling/Model%20Equasions/TripleInvertedPendulum1.txt)
  - [For Matlab](Modeling/Model%20Equasions/TripleInvertedPendulum2.txt)


## Simulation
Explain the simulation processes, tools, and techniques you used. Provide context on what the simulations aim to achieve.

Link to the relevant files and documentation for the simulation section.
- [Simulation Code](path/to/simulation_code)
- [Simulation Documentation](path/to/simulation_documentation)

Include an image or animation that illustrates the simulation results or process .
![Simulation Animation](path/to/simulation_animation.gif)

## Experimental Testing
Describe the experimental testing procedures, including the setup, methodology, and equipment used. Highlight key results and findings.

Link to the relevant files and documentation for the experimental testing section.
- [Experimental Testing Code](path/to/experimental_testing_code)
- [Experimental Testing Documentation](path/to/experimental_testing_documentation)

Include an image or animation related to the experimental testing.
![Experimental Testing Image](path/to/experimental_testing_image.jpg)
