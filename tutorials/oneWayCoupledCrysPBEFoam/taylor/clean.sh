#!/bin/bash
rm -rf [1-9]* \
       processor* \
       postProcessing \
       logs \
       dynamicCode \
       dynamicMesh \
       uniform/time \
       log.*

# 如果你想连网格一起清理，取消下一行注释
# rm -rf constant/polyMesh
