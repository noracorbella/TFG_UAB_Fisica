We tracked the energy of a cellsorting simulation along MCS. 
The noticeable initial increase in energy appears counterintuitive for a Monte Carlo simulation, where energy should decrease towards a minimum.
This initial rise is attributed to the initially unrealistic, square shapes of the cells. 
During the first few Monte Carlo steps (MCS), the volume and surface constraint terms in the Hamiltonian drive the cells to adopt more rounded and realistic shapes.
The total energy plot does not stabilize at a clear minimum, as should be expected. 
This is because CompuCell3D simulations incorporate other factors beyond simple energy minimisation,
like cell motility and adhesion. The system seeks equilibrium reflecting a balance of these forces, rather than a static, minimal energy state.
Even when cell sorting appears complete, cells may continue to exhibit minor rearrangements due to thermal fluctuations.
However, the contact energy does reach a minimum. 
This indicates that the cell sorting process effectively minimises the overall contact energy of the system.
