# connect to jupyter notebook server from xpacc-serv-01

ssh -N -f -L localhost:8888:localhost:8889 $XPACC1
