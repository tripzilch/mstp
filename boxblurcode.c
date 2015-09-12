
    for (int x = 0; x < dimx; ++x) {
      double value = 0.0;
      for (int y = 0; y <= dimy + 2*scale; ++y) {
        value += source(x, y % dimy);
        if (y >= 2*scale) {
          target(x, (y - scale) % dimy) = value / (2*scale+1);
          value -= source(x, (y - 2*scale) % dimy);
        } 
      }
    }
    
    for (int y = 0; y < dimy; ++y) {
      double value = 0.0;
      for (int x = 0; x <= dimx + 2*scale; ++x) {
        value += source(x % dimx, y);
        if (x >= 2*scale) {
          target((x - scale) % dimx, y) = value / (2*scale+1);
          value -= source((x - 2*scale) % dimx, y);
        } 
      }
    }
    
