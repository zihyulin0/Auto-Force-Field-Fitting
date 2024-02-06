# TAFFI - Topology Automated Force Field Interactions
### Welcome to Lin's repo! 
### Unfortunately, all the work I've done as a software developer (not as a researcher) is all under copyright and cannot be shared publicly. 
### This repo is my *ongoing* effort (a lot more unit tests to be added!) to refactor the research code I wrote to have a better data structure and for me and my PhD advisor to release this in the future for free.
- ### For example driver script to see my coding style, see: driver_example.py
- ### For the theory behind this force field framework that we developed, see: https://doi.org/10.1021/acs.jcim.1c00491
- ### If you use any part of the code, please be courteous and cite the paper above!
- ### The data structure hierarchy:
```mermaid
classDiagram
  StructureBase <|--AdjacencyMatrix
  StructureBase <|--AtomTypeBaseOperation
  AdjacencyMatrix<--AtomTypeBaseOperation
  AtomTypeBaseOperation<|--AtomTypes
  AtomTypes<|--Molecule
  class StructureBase{
      +elements
      +geometry
      +q_tot
      +parse_data(data): call static parser
      +parse_xyz(xyz): call static parser
  }
  class AdjacencyMatrix{
      +adj_mat
      +Table_generator()
      +graph_seperation()
      +Dijkstra()
  }
  class AtomTypeBaseOperation{
    +AdjacencyMatrix AdjMat
    +adj_mat: @property
    +graph_seps: @property
    +setter function that sets across AdjMat
    +is_nitro()
    +ring_atom()
    +atom_hash()
    +rec_sum(): ...etc
  }
  class AtomTypes{
    +bond_depth
    +atom_types
    +id_types()
    +find_lewis()
  }
  class Molecule{
    +formal_charge
    +bond_mat
    +atom_types
    +bonds
    +angles
    +dihedrals
    +bond_types
    +angle_types...etc
    +find_modes()
    +write_xyz()
    +update_atom_types():...etc  
  }

```
