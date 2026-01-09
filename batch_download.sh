#!/bin/bash

for idname in 'b7366a5b-c6ca-4ec8-9255-79ccc8f205ac', 'c483e4a0-285f-498e-993f-d3d4dff558a3'; do
	echo "Y" | rscli download "Desktop/Local_repos/streamtemp/Projects/$idname/" --id $idname
done