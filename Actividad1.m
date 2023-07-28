%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

functionNumber = input("Funcion a optimizar: "); %Numero de funcion a ejecutar

f1 = @(x,y) x.*e.^(-x.^2-y.^2); %Funcion 1

f2 = @(x,y) (x-2).^2 + (y-2).^2; %Funcion 2

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Seleccion de padres
function [n] = Selection(apt)
  totalAptitud = sum(apt); %Suma de las aptitudes
  aptitudSize = numel(apt); %Tama�o de las aptitudes
  
  random = rand(); %Numero aleatorio
  selectionSum = 0; %Actual de la suma
  
  for i=1:aptitudSize
    selectionSum = selectionSum + apt(i)/totalAptitud; %Sumar la aptitud actual
    
    if selectionSum >= random
      n = i; %Aptitud seleccionada
      return;
    endif
    
  end
  
  n = aptitudSize; %Seleccionar la ultima aptitud
end

%Cruzar padres
function [y1,y2] = Cross (x1,x2)
  xSize = numel(x1); %Tama�o del padre
  crossPoint = randi([1,xSize]); %Seleccionar punto de cruce
  
  y1 = x1; %Copiar el padre
  y2 = x2; %Copiar el padre
  
  y1(crossPoint:end) = x2(crossPoint:end); %Combinar atributos desde el punto de cruce de x2
  y2(crossPoint:end) = x1(crossPoint:end); %Combinar atributos desde el punto de cruce de x1
  
end

%Mutar hijos
function [y] = Mutation (s,xl,xu)
  [row,column] = size(s); %Tama�o de los hijos
  mutationProbability = 0.01; %Probabilidad de mutaci�n
  y = s; %Copiar hijos
    
  for i=1:column
    for j=1:row
      if rand() < mutationProbability
        y(j,i) = xl(j)+(xu(j)-xl(j))*rand(); %Mutar coordenada de hijo
      endif
      
    end
    
  end
  
end

%Seleccionar funci�n con algoritmo GA
switch functionNumber
  case 1
    xl = -2; %Parte baja en x
    xu = 2; %Parte alta en x
    yl = -2; %Parte baja en y
    yu = 2; %Parte alta en y
    
    [x,y] = meshgrid(xl:0.2:xu,yl:0.2:yu); %Coordenadas a graficar
    z = f1(x,y); %Valores de la funcion evaluada

    contour(x,y,z) %Graficar
    
    %Inicializacion y aptitud
    xl = [-2, -2]'; %Limites inferiores
    xu = [2, 2]'; %Limites superiores
    
    n = 100; %Cantidad de hijos
    d = 2; %Atributos del hijo
    g = 100; %Cantidad de generaciones
    
    x = zeros(d,n); %Inicializar padres en 0
    aptitud = zeros(1,n); %Inicializar aptitudes en 0
    
    %Inicializar padres aleatoriamente
    for i=1:n
      x(:,i) = xl+(xu-xl).*rand(d,1);
    end
    
    %Total de generaciones
    for k=1:g
      %Calcular la aptitud de cada individuo
      for i=1:n
        functionAptitud = f1(x(1,i),x(2,i));
      
        if functionAptitud >= 0
          aptitud(i) = 1/(1+functionAptitud);
        else
          aptitud(i) = 1 + abs(functionAptitud);
        endif
      end
      
      %Conjunto vac�o de hijos
      sons = zeros(d,n);

      %Generar hijos
      for i=1:2:n
        z1 = Selection(aptitud);
        z2 = Selection(aptitud);
        
        while z1 == z2
          z2 = Selection(aptitud);
        end
        
        [son1,son2] = Cross(x(:,z1),x(:,z2)); %Genear dos hijos
        
        sons(:,i) = son1;
        sons(:,i+1) = son2;
        
      end

      sons = Mutation(sons,xl,xu); %Mutar hijos
      
      x = sons; %Los hijos son los padres
            
    end
    
    %Calcular la aptitud de cada individuo
    for i=1:n
      functionAptitud = f1(x(1,i),x(2,i));
    
      if functionAptitud >= 0
        aptitud(i) = 1/(1+functionAptitud);
      else
        aptitud(i) = 1 + abs(functionAptitud);
      endif
    end
    
    %Graficar minimo
    [max,position] = max(aptitud);
    plot(x(1,position),x(2,position),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,position));
    disp(x(2,position));
    disp(f1(x(1,position),x(2,position)));

  case 2
    xl = -5; %Parte baja en x
    xu = 5; %Parte alta en x
    yl = -5; %Parte baja en y
    yu = 5; %Parte alta en y
    
    [x,y] = meshgrid(xl:0.2:xu,yl:0.2:yu); %Coordenadas a graficar
    z = f2(x,y); %Valores de la funcion evaluada

    contour(x,y,z) %Graficar
    
    %Inicializacion y aptitud
    xl = [-5, -5]'; %Limites inferiores
    xu = [5, 5]'; %Limites superiores
    
    n = 100; %Cantidad de hijos
    d = 2; %Atributos del hijo
    g = 100; %Cantidad de generaciones
    
    x = zeros(d,n); %Inicializar padres en 0
    aptitud = zeros(1,n); %Inicializar aptitudes en 0
    
    %Inicializar padres aleatoriamente
    for i=1:n
      x(:,i) = xl+(xu-xl).*rand(d,1);
    end
    
    %Total de generaciones
    for k=1:g
      %Calcular la aptitud de cada individuo
      for i=1:n
        functionAptitud = f2(x(1,i),x(2,i));
      
        if functionAptitud >= 0
          aptitud(i) = 1/(1+functionAptitud);
        else
          aptitud(i) = 1 + abs(functionAptitud);
        endif
      end
      
      %Conjunto vac�o de hijos
      sons = zeros(d,n);

      %Generar hijos
      for i=1:2:n
        z1 = Selection(aptitud);
        z2 = Selection(aptitud);
        
        while z1 == z2
          z2 = Selection(aptitud);
        end
        
        [son1,son2] = Cross(x(:,z1),x(:,z2)); %Genear dos hijos
        
        sons(:,i) = son1;
        sons(:,i+1) = son2;
        
      end
      
      sons = Mutation(sons,xl,xu); %Mutar hijos
      
      x = sons; %Los hijos son los padres

    end
    
    %Calcular la aptitud de cada individuo
    for i=1:n
      functionAptitud = f2(x(1,i),x(2,i));
    
      if functionAptitud >= 0
        aptitud(i) = 1/(1+functionAptitud);
      else
        aptitud(i) = 1 + abs(functionAptitud);
      endif
    end
    
    %Graficar minimo
    [max,position] = max(aptitud);
    plot(x(1,position),x(2,position),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,position));
    disp(x(2,position));
    disp(f2(x(1,position),x(2,position)));

  otherwise
    disp("Numero no valido")

endswitch