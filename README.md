# ocean-cfd

Co-op project at the Bedford Institute of Oceanography (Fisheries and Oceans Canada), looking at how offshore structures stir up the water around them.

There are two things going on in here:

**Wave flume simulations** in OpenFOAM / waves2Foam. Regular waves, a cylinder sitting in the flume, and combined wave plus current cases.

**Passive tracer transport.** This one takes archived LES velocity from a monopile in a steady current and uses it to push a passive tracer around, so you can watch how the wake behind the structure redistributes material through the water column.

## Some results

The cleanest way to see what the cylinder actually does is to subtract the empty flume from the cylinder case. Whatever is left over is purely the structure's doing:

![Regular waves, cylinder vs empty flume](flume_tasks/results/task1_vs_task2_difference.gif)

Same idea for a steady current instead of waves:

![Steady current, cylinder vs no cylinder](flume_tasks/results/task3_vs_task4_difference.gif)

Tracer animations are in `monopile_tracer/*/results/`, though like the flume output they are not tracked here.

## What is where

| | |
|---|---|
| `flume_tasks/` | the four waves2Foam cases and the difference animations |
| `monopile_tracer/` | tracer post-processing of the archived LES |
| `parameter_study/` | earlier wave sweep on a monopile, with its own MANIFEST |
| `archive/` | superseded work, kept for reference |
| `reference/` | third-party cases consulted while learning, not tracked |

Each project directory has its own ReadMe with the detail.

Mesh files, `processor*/`, and the numeric time directories are all gitignored since they are big and regenerable. Run `blockMesh` and `decomposePar` to get them back.

## Built with

* [OpenFOAM](https://www.openfoam.com/) v2412
* [waves2Foam](https://github.com/ogoe/waves2Foam) for wave generation and relaxation zones. Note this is the fork, not upstream, built from `a8d38fd`.
* [olaFlow](https://github.com/phicau/olaFlow) for the earlier monopile work, built from `a35bc3d`.

Python side is just `numpy`, `scipy`, `matplotlib`, and `pandas`.

## A couple of notes

The binned velocity fields and the rendered animations are not tracked here. They run to tens of GB and can be regenerated from the case output.

The tracer starts from a real nitrate profile, Atlantic Zone Monitoring Program station HL5 on the Halifax Line. That data is publicly available from DFO, so it is not redistributed here.
