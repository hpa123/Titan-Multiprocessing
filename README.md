# Titan-Multiprocessing

A computational physics project which simulates the deposition of ions (particles with net electric charge) into the atmosphere of Titan, Saturn's largest moon. Ions in Saturn's magnetospheric plasma collide with Titan's atmosphere as the moon passes through the plasma. The ions deposit energy and induce ionization in the atmosphere. The goal of the project was to demonstrate whether these impacting ions alone could explain observations of free electron density and temperature fluctuations in Titan's atmosphere, which were recorded by the Cassini-Huygens space probe.

*Note: if you try to download and run this code, it will not work, as I have excluded the very large collections of electromagnetic field values which inform the calculations done to track particle trajectories outside of Titan's atmosphere.*

The first important part of the code is the **titanUtil.py** file. This contains the *Particle* class, which is the fundamental unit of the simulation. The instances of this class are *superparticles*, which have statistical weights equivalent to a large quantity of physical particles which enter Titan's atmosphere. Going forward, we will use the terms superparticle and particle interchangeably.

The methods of this class are primarly concerned with calculating relevant physical quantities and storing histories of the (super)particles' trajectories and energies. We also calculate the particles' collision cross-section (essentially, an individual particle's likelihood to collide with atmospheric nitrogen, which changes as a function of the particle's trajectory and energy).

The other two important parts are **mag_traject.py** and **atmos_traject.py**. These deal with particle behaviors outside and inside Titan's atmosphere, respectively. Outside the atmosphere, electromagnetic forces dominate, and the particles are guided by the Lorentz force. Particles can either escape the bounds of the simulation (a box defined under Titan and Grid Attributes in **particle_trajectory_and_energy_dep.py**) or collide with Titan's atmopshere. Those that collide are stored in a seperate file called atmosParticles.npy, which is passed to atmos_traject.py. Particles enter the atmosphere with some initial velocity determined by their trip through the magnetosphere plus the bulk velocity of the plasma (from the particle's perspective, Titan is moving toward it along its orbital trajectory). After that, their fates are determined by their collisions with atmospheric nitrogen.

When all is said and done, we plot energy deposition and ionization as a function of altitude, along with histograms of incident energies and angles.

Both **mag_traject.py** and **atmos_traject.py** contain worker functions which are used with the multiprocessing module to calculate many particle trajectories simultaneously. As such, individual superparticles do not interact in the simulation, as they all exist in their own parallel (processing) universes.
