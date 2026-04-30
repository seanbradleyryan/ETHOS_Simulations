# Context

## Background

You are an expert at RayStation and RayStation scripting. Imported to RayStation are the two CBCT files from `pipeline_setup.m`, the exploded plans, the RTSTRUCT files for the CBCT's, and the frame-of-reference registration between the two images. 

## Goal

Compute and export a dose file for each beam in each plan using a RayStation python script. 

# Technical Details

- The programming language is Python. The script is called `calc_beam_plan_doses.py`. 
- The Linac is Halcyon SN1293. It uses 6FFF photons, has X and Y jaws, and is a dual leaf MLC. 