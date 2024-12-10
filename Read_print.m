function [parametros_S,z0] = Read_print

[file, path] = uigetfile('*.s*p', 'Selecciona un archivo Touchstone (.sNp)');
if isequal(file, 0)
    disp('No se seleccionó ningún archivo. Saliendo del programa.');
    return;
end

% Leer el archivo
filename = fullfile(path, file);
file_readed = fopen(filename, 'r');
if file_readed == -1
    error('No se pudo abrir el archivo.');
end

parametros_S=struct('Impedance',0,'NumPorts',0,'Parameters',0,'Frequencies',0);
n_matriz=0;

while true

    line = fgetl(file_readed);
    

    if line == -1      %final documento
        break;
    end

%% Interpretación de información
    if startsWith(line, '#')
        info = strsplit(line);
        if matches(upper(info{2}),"THZ") == 1
            Hz=1E12;
        end
        if matches(upper(info{2}),"GHZ") == 1
            Hz=1E9;
        end
        if matches(upper(info{2}),"MHZ") == 1
            Hz=1E6;
        end
        if matches(upper(info{2}),"KHZ") == 1
            Hz=1E3;
        end
        if matches(upper(info{2}),"HZ") == 1
            Hz=1;
        end
        z0=str2double(info{6});
    end
    
  if ~startsWith(line, '!') && ~startsWith(line, '#')
        param=1;
        parametros = str2double(strsplit(line));
        matriz_size = sqrt((length(parametros)-2)/2);
        n_matriz=n_matriz+1;
        frequency(n_matriz) = parametros(2)*Hz;
        Matriz_prov = zeros(matriz_size,matriz_size);
        Matriz_S(:,:,n_matriz) = Matriz_prov;

        if matches(upper(info{4}),"RI")
        for i = 1 : 1 : matriz_size
                for k = 1 : 1 : matriz_size
                    param=param+2;
                    Matriz_S(k,i,n_matriz) = parametros(param)+1i*parametros(param+1);
                end
            end
        else
            for i = 1 : 1 : matriz_size
                for j = 1 : 1 : matriz_size
                    param=param+2;
                    RI=parametros(param)*exp(1i*parametros(param+1));
                    Matriz_S(j,i,n_matriz) = RI;
                end
            end
        end

  end

    
end
        parametros_S.Impedance = z0;
        parametros_S.NumPorts = matriz_size;
        parametros_S.Frequencies = frequency;
        parametros_S.Parameters = Matriz_S;

end
