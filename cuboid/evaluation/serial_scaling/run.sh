for x in {1..5}; do
  for y in {1..5}; do
    for z in {1..5}; do
      for f in {1..5}; do
        echo "x=$x, y=$y, z=$z, f=$f"
        $OPENCMISS_REL_DIR/laplace_fortran $OPENCMISS_INPUT_DIR $x $y $z $f
      done
    done
  done
done

