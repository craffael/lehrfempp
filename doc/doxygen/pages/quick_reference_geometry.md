# Quick Reference - Geometry {#quick_reference_geometry}

[TOC]

LehrFEM++ provides an Interface for geometry data.

## Geometry Interface

To get the geometry of an entity

```cpp
for (const lf::mesh::Entity* ent : mesh.Entities(codim)) {

}
```

it can of course also be
