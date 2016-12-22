# strong scaling
while true; do
  for p in 1 2 4 8 16 32 40 48; do
    echo $p
    #mpirun -n $p $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 8 8 8 1 2>&1 | tee -a out${p}.txt
  done
done
