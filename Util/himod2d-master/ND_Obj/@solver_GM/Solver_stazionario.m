function []=Solver_stazionario(obj)
% ATTENZIONE SPERIMENTALE!!!
% L'idea è di fare in modo che funzioni anche in stazionario senza modificare i file che lo facevano andare in tempo dipendente
% Praticamente finito, è necessario creare dei nuovi file.edp senza theta metodo da chiamare al posto dei soliti
if(~obj.assembled)
    display('assembling the problem')
    obj.Assemble;
end
if(obj.solved)
    display('Already solved, use solver_GM.Restart, to compute a new simulation');
else
    obj.solved=true;
    clear raccoltaerr raccoltaerr1 raccoltaerr2 raccoltaerr3 raccoltaerr4
    if(obj.full2d)
        !FreeFem++ -ne -nw ND_full2d_staz.edp
    else
        !FreeFem++ -ne -nw init_staz.edp
        %Inizializzazione variabili ciclo domain decomposition
        errtest=1;
        obj.n_iter(1)=0;
        kap=1;
        
        % ciclo domain decomposition
        while(errtest>obj.toll && obj.n_iter(1)<obj.maxiter)
            %lancio di FreeFem++ per il problema 2D
            !FreeFem++ -ne -nw ND2d_staz.edp
            
            obj.n_iter(1)=obj.n_iter(1)+1;
            file=fopen('testedato2d','r');
            temporaneo=fscanf(file,'%f');
            fclose(file);
            errtest2d=temporaneo(1);
            val_diric_interfaccia=temporaneo(2);
            errtest2dDominio=temporaneo(3);
            if(obj.krpo)
                global w
                global mediaestremi
                w=temporaneo(4);
                mediaestremi=temporaneo(5);
            end
            
            %rilassamento
            %Viene usato da ND1d
            val_diric_interfaccia=obj.thetaDD*val_diric_interfaccia+(1-obj.thetaDD)*obj.u1dkmeno1(1);
            
            %risolvo problema 1D
            ND1d;
            obj.ExportNeum;
            
            % calcolo grandezza test
            errtest=errtest1d*(obj.front1d)+errtest1dDominio*(obj.dom1d)+errtest2d*(obj.front2d)+errtest2dDominio*(obj.dom2d);
            raccoltaerr{1,kap}=errtest;
            raccoltaerr1{1,kap}=errtest1d;
            raccoltaerr2{1,kap}=errtest1dDominio;
            raccoltaerr3{1,kap}=errtest2d;
            raccoltaerr4{1,kap}=errtest2dDominio;
            kap=kap+1;
        end
    end % fine ciclo Domain Decomposition
    
    if(obj.n_iter(1)==obj.maxiter)
        display('Warning: not converging')
        pause;
    end
    obj.SaveTimeStep(1);
end
end
