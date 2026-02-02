// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __PDETOOLS_MESH_Extended_PDE_Mesh_Utilities_HPP
#define __PDETOOLS_MESH_Extended_PDE_Mesh_Utilities_HPP

#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Mesh
{
namespace PDE_Mesh_Utilities
{
  class PDE_Domain_2D_Collection
  {
    public:
      std::vector<Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D> domains_2D;
      std::vector<Eigen::Matrix3d> domains_rotation;
      std::vector<Eigen::Vector3d> domains_translation;
  };

  class Extended_MeshGeometricData2D
  {
    public:
      Gedim::MeshUtilities::MeshGeometricData2D mesh_geometric_data;
      std::vector<unsigned int> cell2Ds_domain_position;
  };

} // namespace PDE_Mesh_Utilities
} // namespace Mesh
} // namespace PDETools
} // namespace Polydim

#endif
