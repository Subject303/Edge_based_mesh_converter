#!/bin/bash

# case_file="half_domain_slender.vtk"
# case_file="full_domain.vtk"
case_file="slender-no-wake-cour.vtk"
# case_file="slender-no-wake-acc.vtk"
# case_file="star1.case"

ls

. "./run_edge_based_mesh_converter.sh" $case_file
