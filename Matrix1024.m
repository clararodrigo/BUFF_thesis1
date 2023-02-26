classdef Matrix1024 < BaseTransducer
    
    properties (Access = public)
        full_aperture_tx;
        full_aperture_rx;
        SpAp;
        f0 = 7.8125*MHz;
    end
    
    methods
        function self = Matrix1024()
            self.name = "Matrix1024";
            self.description = "Matrix Array transducer";

            self.full_aperture_tx = GridAperture( ...
                [32, 35],            ... % Count
                0.275*mm,           ... % Width
                0.275*mm,           ... % Height
                [0.025*mm, 0.025*mm],   ... % Kerf % should be 0.025
                1,                  ... % Nx
                1                   ... % Ny
            );

            % remove gaps
            self.full_aperture_tx.elements = self.full_aperture_tx.elements([1:8,10:17,19:26,28:35],: );
            for ei = 1:self.full_aperture_tx.n_elements
                self.full_aperture_tx.elements(ei).number = ei;
            end
            self.full_aperture_tx.build() % rebuild
            
            % Create RX aperture as copy of TX aperture
            self.full_aperture_rx = self.full_aperture_tx.deep_copy();
            self.full_aperture_rx.fII_id = [];
            self.full_aperture_rx.build() % rebuild
            
            % Set Impulse responses and excitations
            self.full_aperture_tx.impulse_response = F0BWImpulseResponse([2.4290*MHz 4.5110*MHz], self.f0, 100);
            self.full_aperture_tx.apply_waveforms();
            self.full_aperture_tx.bw = [2.4290*MHz 4.5110*MHz];
            self.full_aperture_tx.excitation.f0 = self.f0;
            
            self.full_aperture_rx.impulse_response = F0BWImpulseResponse([2.4290*MHz 4.5110*MHz], self.f0, 100);
            self.full_aperture_rx.apply_waveforms();
            self.full_aperture_rx.bw = [2.4290*MHz 4.5110*MHz];
            self.full_aperture_rx.excitation.f0 = self.f0;
                        
            self.tx_aperture = self.full_aperture_tx;
            self.rx_aperture = self.full_aperture_rx;
        end
        
        function initial_load_apodizations(self, fname)
            SpAp = load(fname, 'CompRndApod').CompRndApod; 
            self.SpAp = SpAp;
        end
        
        function set_receive_apodization(self, set_n, apo_n)
            self.rx_aperture.elem_apodization = self.SpAp(set_n, apo_n, :);
            self.rx_aperture.apply_apodization();
        end
        
        function set_transmit_apodization(self, set_n)
            self.tx_aperture.elem_apodization = self.SpAp(set_n, 1, :);
            self.rx_aperture.apply_apodization();
        end
    end
end      
