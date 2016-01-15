% en los resultados, el nombre es de la forma:
%  In_#_com  donde # es el indice de la incidencia y 
%                  com es el nombre del componente

% se carga una variable de nombre de la forma:
%%%  inci_a_b  donde   a   es el numero de incidencia
%%%                    b   es un indice del componente
% El nombre de la variable es 'var'

% cada renglon es un receptor
% cada columna es un dato en el tiempo

%% 
In_2_szz
var = real(var);
%
maxAmp = max(max(var));
offset = maxAmp*1.5;
[nPts,nmuestras] = size(var);

figure; hold on
for ipt = 1:4
    plot(var(ipt,:)+(ipt-1)*offset,'k')
end

%%
figure; hold on
for ipt = 5:nPts
    plot(var(ipt,:)+(ipt-1)*offset,'k')
end 
