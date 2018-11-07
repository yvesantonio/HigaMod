function coeff_output = interp_interface( coeff_input, ril_input, ril_output, mb_input, mb_output,pesi,bc)
%
%   function coeff_output = interp_interface( coeff_input, ril_input, ril_output, mb_input, mb_output, pesi,bc )
%
%	riceve coeff_input: un vettore di coefficienti lungo come la base modale in input ossia la base modale che si affaccia alla frontiera
%          ril_input  : il rilevamento va passanto gi?? valutato nei nodi di quadratura

% TODO aggiungere ottimizzazione se le basi sono uguali c'?? da proiettare solo il rilevamento, se anche quello ?? uguale niente da fare.
%      probabilmente conviene farlo in un altro punto del codice prima di chiamare coeff_output una roba del tipo se serve chiama interp_interface
%      altrimenti non lo chiamare


m_out = size( mb_output , 2 ); % numero elementi in base dominio di interesse (sx)
m_in  = size( mb_input  , 2 ); % numero elementi in base dominio che si affaccia (dx)
coeff_output=zeros(m_out, 1 );

if (strcmp(bc,'dir')||strcmp(bc,'rob'))
    ril   = ril_input - ril_output;  % differenza tra i rilevamenti:
end

for i = 1 : m_out
    
    if (strcmp(bc,'dir') || strcmp(bc,'rob'))
        coeff_output(i) = integrate( (ril.*mb_output(:,i))', pesi ); %ATTENZIONE manca jacobiano (L)
    end

    for j = 1 : m_in
        coeff_output(i) = coeff_output(i) + coeff_input(j)*integrate((mb_output(:,i).*mb_input(:,j))',pesi);%ATTENZIONE manca jacobiano (L)
    end    
end
