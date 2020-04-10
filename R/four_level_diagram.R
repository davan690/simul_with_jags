grViz(" 
  digraph CFA {
    # Multiple level nodes
    node [shape = ellipse, color=CornflowerBlue]
    a [label = '@@1']; 
    b [label = '@@2']; 
    c [label = '@@3']; 
    d [label = '@@4'];
    {rank = same; b; d}

    # Terminal branch nodes
    node [shape = box, color = Crimson] 
    e [label = 'Model 2'];
    f [label = 'Model 4'];
    g [label = 'Model 5'];
    h [label = 'Model 7'];
    i [label = 'Model 8'];
    {rank = same; e; f; g; h; i}

    # Connect nodes with edges and labels
    a -> b [label = 'Condition 1a']
    a -> d [label = 'Condition 1b'] 
    b -> e [label = 'Condition 2a'] 
    b -> c [label = 'Condition 2b']
    c -> f [label = 'Condition 3a']
    c -> g [label = 'Condition 3b']
    d -> h [label = 'Conddition 4a'] 
    d -> i [label = 'Condition 4b'] 
  }

[1]: 'Split 1' 
[2]: paste0('Model 1\\n Split 2') 
[3]: paste0('Model 3\\n Split 3') 
[4]: paste0('Model 6\\n Split 4') 
")