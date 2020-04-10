library("DiagrammeR")
grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = circle,
        fontname = Helvetica]
  rats
  
    # several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  seed

  # several 'edge' statements
  rats->seed; rats->rats; seed->seed; seed->rats
}
")
