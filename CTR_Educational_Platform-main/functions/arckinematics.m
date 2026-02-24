function T = arckinematics(k, phi, l)
            if k > 0
                R = [cos(phi)*cos(k*l), -sin(phi), cos(phi)*sin(k*l);
                    sin(phi)*cos(k*l), cos(phi), sin(phi)*sin(k*l);
                    -sin(k*l), 0, cos(k*l)];
                
                p = [(cos(phi) * (1-cos(k*l))) / k;
                            (sin(phi) * (1-cos(k*l))) / k;
                            sin(k*l) / k];
            else 
                R = [cos(phi), -sin(phi), 0;
                    sin(phi), cos(phi), 0;
                    0 0 1];
                p = [0; 0; l;];
                
            end 
            T = [R p;
                0 0 0 1];

        end