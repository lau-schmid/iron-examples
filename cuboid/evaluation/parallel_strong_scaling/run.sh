# strong scaling


for p in 32 16 8 4 2 1; do
  echo $p
  mpirun -n $p $OPENCMISS_REL_DIR/cuboid $OPENCMISS_INPUT_DIR 8 8 8 1 2>&1 | tee -a out${p}.txt
done
