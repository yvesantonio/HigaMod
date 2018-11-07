function []=Solve(obj)
% function []=Solve(obj)
%
% Risolve il problema, se questo non è stato assemblato lo assembla
% Se è già stato risolto occorre chiamare solver_GM.Restart
% L'output è accessibile tramite le funzioni di Plot 
% e le variabili Solution e n_iter.
% Solution è una struttura dati in cui il primo indice corrisponde
% all'istante temporale (1 equivale al tempo iniziale)
% il secondo rappresenta la zona 1=2d, 2=1d.
% mette solved a vero
if(obj.timedipendent)
if(~obj.assembled)
    display('assembling the problem')
    obj.Assemble;
end
if(obj.solved)
    display('Already solved, use solver_GM.Restart, to compute a new simulation');
else
    obj.solved=true;
for i=1:obj.n_passi
    if(~obj.first)
        clear raccoltaerr raccoltaerr1 raccoltaerr2 raccoltaerr3 raccoltaerr4
    end
    if(obj.full2d)
            !FreeFem++ -ne -nw ND_full2d.edp
        
        if(~obj.first)
            obj.t=obj.t+obj.dt;
            obj.SaveTimeStep(i);
        else
            obj.first=0;
            obj.ExportFirst;
            obj.SaveTimeStep(i);
        end
    else
        %Inizializzazione variabili ciclo domain decomposition
        errtest=1;
        obj.n_iter(i)=0;
        kap=1;
        
        % ciclo domain decomposition
        while(errtest>obj.toll && obj.n_iter(i)<obj.maxiter)
            %lancio di FreeFem++ per il problema 2D
                !FreeFem++ -ne -nw ND2d.edp
            
            if(obj.first)%se non � la prima volta plotto il dato iniziale
                errtest=obj.toll/2;
                obj.u1d=zeros(length(obj.u1dold)*2,1);
                for h=1:length(obj.u1dold)
                    obj.u1d(2*h-1)=obj.u1dold(h);
                    obj.u1d(2*h)=obj.u1dold(h);
                end
                %inizializzazione loop spaziale prima volta.
                !cp u2dold u2dkmeno1
            else
                obj.n_iter(i)=obj.n_iter(i)+1;
                [file,messaggio]=fopen('testedato2d','r');
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
                % viene usato da ND1d
                val_diric_interfaccia=obj.thetaDD*val_diric_interfaccia+(1-obj.thetaDD)*obj.u1dkmeno1(1);
                
                %risolvo problema 1D
                ND1d;
                
                obj.ExportNeum;
                
                % calcolo grandezza test
                errtest=errtest1d*(obj.front1d)+errtest1dDominio*(obj.dom1d)+errtest2d*(obj.front2d)+errtest2dDominio*(obj.dom2d);
                raccoltaerr{i,kap}=errtest;
                raccoltaerr1{i,kap}=errtest1d;
                raccoltaerr2{i,kap}=errtest1dDominio;
                raccoltaerr3{i,kap}=errtest2d;
                raccoltaerr4{i,kap}=errtest2dDominio;
                kap=kap+1;
            end
        end % fine ciclo Domain Decomposition
        
        if(~obj.first)
            %esporto condizioni di bordo problema 2D
            neum_per_2d=fopen('neum','r');
            aus2=fscanf(neum_per_2d,'%f');
            obj.val_neum_old=aus2(1);
            fclose(neum_per_2d);
            obj.ExportNeum;
        end
        
        if(obj.n_iter(i)==obj.maxiter)
            display('Warning: not converging')
            pause;
        end
        
        %aggiornamento soluzione passo precedente
        if(~obj.first)
            obj.t=obj.t+obj.dt;
            %aggiornamento temporale  1d
            obj.u1dold=u1dnew;
            %aggiornamento temporale 2d
            !cp u2dkmeno1 u2dold
            obj.SaveTimeStep(i);
            
        else
            obj.first=0;
            obj.ExportFirst;
            obj.SaveTimeStep(i);
        end
    end
end

end
else
obj.Solver_stazionario;
end
end
