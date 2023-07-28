%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

methodNumber = input("Metodo a utilizar: "); %Metodo a utilizar

load data

m = size(X,1); %Tama�o de datos
g = @(w) w(1) + w(2)*X + w(3)*X.^2 + w(4)*X.^3 + w(5)*X.^4;	%Polinomio
f = @(w) (0.5/m)*sum((Y - g(w)).^2); %Funci�n

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Funcion de recombinaci�n
function y = Recombination (x1,x2)
    y = 0.5*(x1+x2);
end

G = 30; %Numero de generaciones
mu = 150; %Numero de padres
D = 5; %Tama�o de valores
lambda = 180; %Numero de hijos

%Estrategia Evolutiva
  
plot(X,Y,'ro','LineWidth',2,'MarkerSize',10) %Graficar datos

xl = [-8 -3 -3 -1 -0.5]'; %Parte baja
xu = [8 3 3 1 0.5]'; %Parte alta

x = zeros(D,mu+lambda); %Inicializar en ceros los padres
sigma = zeros(D,mu+lambda); %Inicializar en ceros los sigmas
fitness = zeros(1,mu+lambda); %Inicializar en ceros el fitness

%Inicializar aleatoriamente
for i=1:mu
  x(:,i) = xl+(xu-xl).*rand(D,1);
  sigma(:,i) = 0.5*rand(D,1);
  fitness(i) = f(x(:,i));
end

%Generaciones
for k=1:G
  %Crear hijos
  for i=1:lambda
    %Elegir padres
    r1 = randi([1 mu]);
    r2 = randi([1 mu]);
    while r1==r2
        r2 = randi([1 mu]);
    end
    
    %Recombinar padres
    x(:,mu+i) = Recombination(x(:,r1),x(:,r2));
    sigma(:,mu+i) = Recombination(sigma(:,r1),sigma(:,r2));
    
    %Calcular normalizaci�n
    r = normrnd(0,sigma(:,mu+i),[D 1]);
    x(:,mu+i) = x(:,mu+i) + r;
  
    %Evaluar fitness
    fitness(mu+i) = f(x(:,mu+i));
    
  end
  
  switch methodNumber
    %(mu + lamda)-ES
    case 1
      %Ordenar hijos y padres
      [~,I] = sort(fitness);
      x = x(:,I);
      sigma = sigma(:,I);
      fitness = fitness(I);
    
    %(mu,lambda)-ES
    case 2
      %Inicializar hijos
      y = zeros(D,lambda);
      ySigma = zeros(D,lambda);
      yFitness = zeros(1,lambda);
      
      %Obtener hijos
      for j=1:lambda
        y(:,j) = x(:,mu+j);
        ySigma(:,j) = sigma(:,mu+j);
        yFitness(:,j) = fitness(1,mu+j);
      end
      
      %Ordenar hijos
      [~,I] = sort(yFitness);
      y = y(:,I);
      ySigma = ySigma(:,I);
      yFitness = yFitness(:,I);
      
      %Copiar hijos ordenados
      for j=1:lambda
        x(:,j) = y(:,j);
        sigma(:,j) = ySigma(:,j);
        fitness(:,j) = yFitness(:,j);
      end
      
      %Borrar padres
      for j=1:lambda
        x(:,mu+j) = zeros(D,1);
        sigma(:,mu+j) = zeros(D,1);
        fitness(:,mu+j) = zeros(1,1);
      end

  endswitch
  
end

%Graficar minimo
plot(X,g(x(:,1)),'b-','LineWidth',3,'MarkerSize',10)

%Mostrar resultados
disp("Minimos");
disp(x(:,1));
disp("Error");
disp(f(x(:,1)));