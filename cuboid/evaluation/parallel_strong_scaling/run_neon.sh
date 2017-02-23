# strong scaling
  for p in 144 72 48 36 24 18 16 12 9 8 6 4 3 2 1; do
    echo $p
    #mpirun -n $p $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 12 12 12 1 2>&1 | tee -a out${p}.txt
  done
