function skew_w = skew(omega)
    skew_w = [0 -omega(3) omega(2);
              omega(3) 0 -omega(1);
              -omega(2) omega(1) 0];
end 