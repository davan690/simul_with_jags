grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  rats; mice
  
    # several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  seed
  
  # several 'node' statements
  node [shape = circle,
        fontname = Helvetica]
  stoats

  # several 'edge' statements
  rats->mice mice->rats stoats->mice stoats->rats rats->seed mice->seed
}
")