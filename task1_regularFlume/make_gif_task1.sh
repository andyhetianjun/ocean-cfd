#!/bin/bash
dirfold="New Folder"
fprefix="task1_animation"
figs="${dirfold}/${fprefix}*.png"
moviefile="task1_birdseye_animation.gif"
convert -delay 10 -loop 0 ${figs} ${moviefile}
