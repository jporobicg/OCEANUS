# OCEANUS: Ocean Circulation and Transport Analysis Package

OCEANUS is a comprehensive Python package designed for analyzing ocean circulation patterns and computing transport fluxes between spatial boxes in marine environments. The package provides tools for processing hydrodynamic model outputs, calculating inter-box transport, computing state variables, and applying mass balance corrections to ensure physical consistency.

## Transport Between Boxes Using Layers

OCEANUS calculates transport between spatial boxes by integrating velocity fields across predefined depth layers and along face boundaries. The package processes hydrodynamic model outputs (u and v velocity components) and interpolates these fields onto face integration points using linear interpolation methods. For each face connecting adjacent boxes, OCEANUS computes the transport flux by integrating the velocity component normal to the face across the face length and through each depth layer.

The transport calculation follows the fundamental relationship:

$$T_{i,j,k} = \int_{face} \int_{layer} \vec{v} \cdot \hat{n} \, dA \, dz$$

where $T_{i,j,k}$ represents the transport through face $i$ at time $j$ and layer $k$, $\vec{v}$ is the velocity vector, $\hat{n}$ is the unit normal vector to the face, and the integration is performed over the face area and layer depth. OCEANUS handles coordinate transformations, accounts for Earth's curvature through latitude corrections, and manages missing data through robust interpolation schemes.

## Variables Per Box and Layer

OCEANUS computes and retrieves state variables for each spatial box and depth layer by averaging hydrodynamic model outputs over the box volumes. The package processes three-dimensional fields such as temperature, salinity, and vertical velocity, then spatially averages these variables within each box geometry and vertically averages across each predefined depth layer. This approach provides representative values for each box-layer combination while preserving the vertical structure of the ocean.

The spatial averaging process in OCEANUS follows:

$$\bar{\phi}_{i,k} = \frac{1}{V_{i,k}} \int_{V_{i,k}} \phi(x,y,z) \, dV$$

where $\bar{\phi}_{i,k}$ is the averaged variable $\phi$ for box $i$ and layer $k$, $V_{i,k}$ is the volume of the box-layer combination, and the integration is performed over the entire box volume. OCEANUS handles complex box geometries, manages partial box coverage within the hydrodynamic domain, and provides comprehensive diagnostics for data quality assessment.

## Mass Balance Correction

OCEANUS applies a mass balance correction algorithm to ensure that the net transport into and out of each box equals zero, maintaining the fundamental principle of mass conservation. The correction algorithm identifies boxes with non-zero net transport and applies adjustments to faces with missing or zero flux values (NaN or zero transport values). This approach preserves the original transport data while ensuring physical consistency across the entire domain.

The mass balance correction in OCEANUS follows the constraint:

$$\sum_{j \in \text{faces of box } i} T_{j,k} = 0$$

for each box $i$ and layer $k$, where $T_{j,k}$ represents the transport through face $j$ in layer $k$. When this constraint is violated, OCEANUS calculates the imbalance $\Delta T_{i,k} = \sum_{j} T_{j,k}$ and applies a correction $C_{j,k} = -\Delta T_{i,k}$ to an appropriate empty face connected to the imbalanced box. The algorithm prioritizes faces with missing data and ensures that corrections are applied in a physically meaningful manner, maintaining the overall transport structure while achieving mass balance.

---

*For detailed implementation information and usage examples, please refer to the individual module documentation and example scripts provided with the OCEANUS package.*
