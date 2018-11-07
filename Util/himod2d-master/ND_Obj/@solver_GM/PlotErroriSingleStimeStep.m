function PlotErroriSingleStimeStep(obj,num)
%function PlotErroriSingleStimeStep(obj,num)
%
%
%       La funzione stampa ai vari intervalli di tempo gli errori DD
    if((~first))
        figure(2)
        %             subplot(6,1,1);
        bas=1:length(obj.raccoltaerr{num});
        %             semilogy(bas,raccoltaerr);
        %             xlabel('errore totale');

        subplot(2,1,1);
        loglog(bas,obj.raccoltaerr1{num});
        xlabel('errore sul dato di neumann calcolata nel 1d');
        title(['t = ',num2str(t)]);
        xlim([0,maxiter]);
        %              subplot(6,1,3);
        %              semilogy(bas,raccoltaerr2);
        %              xlabel('errore sul dominio 1d');
        %
        subplot(2,1,2);
        loglog(bas,obj.raccoltaerr3{num});
        xlabel('errore 2D sulla frontiera 1d/2d');
        xlim([0,maxiter]);
        %
        %             subplot(6,1,5);
        %             semilogy(bas,raccoltaerr4)
        %             xlabel('errore 2D sul dominio');
    end
end