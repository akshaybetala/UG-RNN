#!/usr/bin/env bash
DATASET=$1


sudo python setup.py install && python -m ugrnn.main --train --model_no=0 --contract_rings
python -m ugrnn.main --train --model_no=1 --contract_rings
python -m ugrnn.main --train --model_no=2 --contract_rings
python -m ugrnn.main --train --model_no=3 --contract_rings
python -m ugrnn.main --train --model_no=4 --contract_rings
python -m ugrnn.main --train --model_no=5 --contract_rings
python -m ugrnn.main --train --model_no=6 --contract_rings
python -m ugrnn.main --train --model_no=7 --contract_rings
python -m ugrnn.main --train --model_no=8 --contract_rings
python -m ugrnn.main --train --model_no=9 --contract_rings
python -m ugrnn.main --train --model_no=10 --contract_rings
python -m ugrnn.main --train --model_no=11 --contract_rings
python -m ugrnn.main --train --model_no=12 --contract_rings
python -m ugrnn.main --train --model_no=13 --contract_rings
python -m ugrnn.main --train --model_no=14 --contract_rings
python -m ugrnn.main --train --model_no=15  --contract_rings
python -m ugrnn.main --train --model_no=16 --contract_rings
python -m ugrnn.main --train --model_no=17 --contract_rings
python -m ugrnn.main --train --model_no=18 --contract_rings
python -m ugrnn.main --train --model_no=19  --contract_rings


