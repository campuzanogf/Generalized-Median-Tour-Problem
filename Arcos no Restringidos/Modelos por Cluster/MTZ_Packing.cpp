/************************************************************************************
* Este es un programa en C++ code para el TSP, formulacion con SECs, que se separan
* mediante el algoritmo de corte mínimo de Stoer y Wagner
************************************************************************************/

#include <iostream> // flujo de entrada y salida
#include <fstream> // flujo de archivos
#include <string> // cadenas de caracteres
#include <cstdlib> // rutinas de C, como rand()
#include <vector> // vectores STL
#include <set> // conjuntos STL
#include <algorithm> // algoritmos STL
#include <cmath> // matematica de C
#include <ilcplex/ilocplex.h> // Librería de CPLEX

#include <stdio.h>
#include <time.h>
#include <ctime>


using namespace std; // Espacio de nombres basico, para que todo funcione
ILOSTLBEGIN // Espacio de nombres de CPL

struct Data // Estructura para almacenar los datos del problema
{
	int cardN;                 // cantidad de nodos
	IloArray<IloNumArray> C;
	IloNumArray SS;// Guardo cluster del nodo casillero
	IloNumArray id;// matriz de costos
	IloNumArray x_cord;
	IloNumArray y_cord;
	double OF;                 // valor (optimo) de la FO
	double time_solve;         // tiempo de resolucion
	char data_filename[255];   // archivo de entrada
	char output_filename[255]; // archivo de salida
	IloArray<IloNumArray> Cluster;         // nodos de cada cluster	
	int clusterN;              // numero de clusteres
};


bool Load_data(IloEnv &env, char data_filename[255], Data &instance); // carga los datos de la instancia
bool Create_MIP(Data &instance, IloModel &Model, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y, IloIntVarArray &t);  // crea el modelo instanciado
bool MIP_solve(IloModel &Model, IloCplex &Cplex, Data &instance, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y, IloIntVarArray &t); // resuelve la instancia
void Guardar(char *nombre, char resultados[50], IloCplex &Cplex, double &time, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y);
void Graficar_yed(char *nombre, Data &instance, IloCplex &Cplex, IloArray<IloBoolVarArray> &X, IloIntVarArray &t, IloArray<IloBoolVarArray> &y, float &ZZ1, float &ZZ2);

// Elementosde trabajo
double eps_rhs; // Tolerancia para lados derechos
int SEC = 0;

class FacilityCallback : public IloCplex::Callback::Function {
private:

	/* Variables representing travel i -> j */
	IloArray<IloBoolVarArray> X;

	/* Variables representing asigment i -> j */
	IloArray<IloBoolVarArray> y;

	/* instancia cargada */
	Data instance;

public:
	/* Constructor with data */
	FacilityCallback(const IloArray<IloBoolVarArray> &_X, const IloArray<IloBoolVarArray> &_y, const Data &_instance) : X(_X), y(_y), instance(_instance) {}

	inline void
		SeparateAlgorithmCuts(const IloCplex::Callback::Context &context) const
	{

		IloArray<IloNumArray> X_Sol(X.getEnv(), instance.cardN);
		for (int i = 0; i < instance.cardN; ++i)
		{
			X_Sol[i] = IloNumArray(X.getEnv(), instance.cardN);
			for (int j = 0; j < instance.cardN; ++j)
			{
				X_Sol[i][j] = context.getCandidateValue(X[i][j]);
			}
		}

		IloArray<IloNumArray> Y_Sol(y.getEnv(), instance.cardN);
		for (int i = 0; i < instance.cardN; ++i)
		{
			Y_Sol[i] = IloNumArray(y.getEnv(), instance.cardN);
			for (int j = 0; j < instance.cardN; ++j)
			{
				Y_Sol[i][j] = context.getCandidateValue(y[i][j]);
			}
		}

		/*
		cout << endl << "Valores de X[i][j]" << endl;
		for (int i = 0; i < instance.cardN; i++)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				if (X_Sol[i][j] > eps_rhs)
					cout << endl << "X[" << i << "][" << j << "] = " << X_Sol[i][j];
			}
		}

		cout << endl << "Valores de Y[i][j]" << endl;
		for (int i = 0; i < instance.cardN; i++)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				if (Y_Sol[i][j] > eps_rhs)
					cout << endl << "Y[" << i << "][" << j << "] = " << Y_Sol[i][j];
			}
		}
		getchar();
		getchar();
		*/

		// Algoritmo de Separacion //
		vector <int> S1; //Cilo identificado
		vector <int> S2; //Nodos por revisar
		vector <int> S3; //Complemento de S1
		vector <int> S4; //Vector de cluster, al cual pertenecen los nodos(identificados)
		vector <int> S5; //Vector de cluster complemento, alcual pertenecen los nodos identificados
		int l = 0; //contados que se utiliza en los cilos de S1, S2 y S3

		//Valor inicial para S2
		for (int i = 0; i < instance.cardN; i++)
			S2.push_back(i);

		int k = 0; //contador de ciclos
		k = 0;
		while (S2.size() != 0)
		{
			//Valores Iniciales
			S1.clear();
			S1.push_back(S2[0]);

			//Identifico ciclo en S1
			for (int i = 0; i < S1.size(); i++)
			{
				for (int j = 0; j < S2.size(); j++)
				{
					l = 0;
					if (S1[i] != S2[j])
					{
						if (X_Sol[S1[i]][S2[j]] > eps_rhs || X_Sol[S2[j]][S1[i]] > eps_rhs || Y_Sol[S2[j]][S1[i]] > eps_rhs) //si encuentro un arco que conecte con nodo S1[i]
						{
							for (int p = 0; p < S1.size(); p++)//me aseguro que el nodo ya no esté en el ciclo
							{
								if (S1[p] == S2[j])
									l++;
							}

							if (l == 0) // si no esta en el ciclo, lo agrego
							{
								S1.push_back(S2[j]);
							}
						}
					}
				}
			}

			//Creo el complemento de S1 en S2
			S3.clear();
			for (int i = 0; i < S2.size(); i++)
			{
				l = 0;
				for (int j = 0; j < S1.size(); j++)
					if (S1[j] == S2[i])
						l++;

				if (l == 0)
					S3.push_back(S2[i]);
			}

			S2.clear();
			for (int i = 0; i < S3.size(); i++)
				S2.push_back(S3[i]);

			//Creo en S3 el complemento de la red para S1
			S3.clear();
			for (int i = 0; i < X.getSize(); i++)
			{
				l = 0;
				for (int j = 0; j < S1.size(); j++)
					if (S1[j] == i)
						l++;
				if (l == 0)
					S3.push_back(i);
			}

			//en S4 guardo los cluster identificados en el subtour
			S4.clear();
			S4.push_back(instance.SS[S1[0]]);
			for (int i = 0; i < S1.size(); i++)
			{
				l = 0;
				for (int j = 0; j < S4.size(); j++)
				{
					if (S4[j] == instance.SS[S1[i]])
						l++;
				}
				if (l == 0)
					S4.push_back(instance.SS[S1[i]]);
			}

			// En S5 creo complemento de S1 en el cluster
			S5.clear();
			for (int i = 0; i < instance.Cluster[S4[0]].getSize(); i++)
			{
				l = 0;
				for (int j = 0; j < S1.size(); j++)
					if (instance.Cluster[S4[0]][i] == S1[j])
						l++;
				if (l == 0)
					S5.push_back(instance.Cluster[S4[0]][i]);
			}

			//elimino la ruta que se identifica el deposito
			if (k == 0)
			{
				S1.clear();
				S3.clear();
				S4.clear();
				S5.clear();
			}

			/*solo para graficar*/
			/*
			 cout << endl << "Ciclo identificado" << endl;
			 for (int i = 0; i < S1.size(); i++)
			 {
			 cout << " " << S1[i]+1;
			 }

			 cout << endl << "clusters identificados" << endl;
			 for (int i = 0; i < S4.size(); i++)
			 {
			 cout << " " << S4[i];
			 }

			 cout << endl << "complemento de clusters identificados" << endl;
			 for (int i = 0; i < S5.size(); i++)
			 {
			 cout << " " << S5[i];
			 }
			 */

			 //Condicion de Entrada: Ciclo cerrado, que no sea de un elemento y que no sea del tamaño de la red (TSP)
			if (S1.size() >= 2)
			{
				/*
				cout << endl << "cluster identificados: " << endl;
				for (int p = 0; p < S4.size(); p++)
				{
					cout << endl;
					for (int j = 0; j < instance.Cluster[S4[p]].getSize(); j++)
					{
						cout << " " << instance.Cluster[S4[p]][j];
					}
				}
				cout << endl;
				*/

				//Corte Packing
				
				IloExpr Edge(X.getEnv());
				IloExpr Assig(y.getEnv());
				IloNum S = 0;
				for (int i = 0; i < S1.size(); i++)
				{
					S++;
					for (int j = 0; j < S1.size(); j++)
					{
						if (S1[i] != S1[j])
						{
							Edge += X[S1[i]][S1[j]];
							Assig += y[S1[i]][S1[j]];
						}
					}
				}
				context.rejectCandidate(Edge + Assig <= S - 1);
				//cout << endl << "Corte: " << endl;
				//cout << endl << Edge Assig << " <= " << S - 1 << endl;
				Edge.end();
				Assig.end();
				

				SEC++;
			}

			//ajusto criterio de salida
			if (S1.size() == instance.cardN)
				S2.clear();

			//Aumento contador
			k++;
			// cout << endl << "Al final del While" << endl;
		}//while(S2.size() != 0) mientras no se hayan revisado todos los nodos de la red

		//Cierro Vector con valores
		for (int i = 0; i < instance.cardN; ++i)
		{
			X_Sol[i].end();
		}
		X_Sol.end();

		for (int i = 0; i < instance.cardN; ++i)
		{
			Y_Sol[i].end();
		}
		Y_Sol.end();

	}

	// This is the function that we have to implement and that CPLEX will call
	// during the solution process at the places that we asked for.
	virtual void invoke(const IloCplex::Callback::Context &context);

	/// Destructor
	virtual ~FacilityCallback();
};

/* Implementation of the invoke function */
void
FacilityCallback::invoke(const IloCplex::Callback::Context &context)
{
	if (context.inCandidate())
	{
		SeparateAlgorithmCuts(context);
	}
}

/// Destructor
FacilityCallback::~FacilityCallback()
{

}



int main(int argc, char **argv)
{
	IloEnv env; // Crear el entorno
	try // ejecutar lo siguiente, alegar si no es el caso
	{
		if (argc == 3) // si se pasan los comandos de la forma 'ejecutable archivo_entrada archivo_salida'
		{
			Data instance;	// crear donde se almacenan los datos

			if (Load_data(env, argv[1], instance)) // probar a cargar los datos
			{
				IloModel Model(env); // crear modelo vacio asociado al entorno

				IloCplex Cplex(Model); // crear ejecucion de CPLEX de este modelo

				/* Crear variables de decision */
				IloArray<IloBoolVarArray> X(env, instance.cardN); // X[i][j]
				IloArray<IloBoolVarArray> y(env, instance.cardN); // y[i][j]
				IloIntVarArray t(env, instance.clusterN); // t[i]

				if (Create_MIP(instance, Model, X, y, t)) // crear el modelo de la instancia
				{
					if (MIP_solve(Model, Cplex, instance, X, y, t)) // resolverlo
					{
						cout << "Estatus en palabras:" << Cplex.getCplexStatus() << endl;

						IloInt inf = Cplex.getCplexStatus();
						cout << "Estatus en numeros: " << inf << endl;

						if (Cplex.getCplexStatus() == CPX_STAT_UNBOUNDED)
						{
							cout << "Problema no acotado" << endl;
						}
						else if (Cplex.getCplexStatus() == CPX_STAT_INFEASIBLE)
						{
							cout << "Problema infactible" << endl;
						}
						else // si la solucion es factible, e idealmente optima
						{

							string instancia; //guarda nombre de la instancia
							int i = 0;
							while (argv[1][i] != '.')
							{
								instancia.push_back(argv[1][i]);
								i++;
							}

							double time = Cplex.getTime();

							//Rsulevo relajacion lineal
							instance.OF = Cplex.getObjValue(); // guardar valor objetivo
							cout << "Instancia: " << instancia << endl; // que se calcula y muestra
							cout << "Tiempo de procesamiento (CPU Time): " << Cplex.getTime() << endl; // que se calcula y muestra
							cout << "Tiempo de procesamiento (wall clock time): " << time << endl; // que se calcula y muestra
							cout << "Solucion optima del problema Z = " << Cplex.getObjValue() << endl;

							//guardar valores
							Guardar(argv[1], argv[2], Cplex, time, X, y);

							//creo valores de funciones objetvos
							float ZZ1;
							ZZ1 = 0;
							for (int i = 0; i < instance.cardN; ++i)
								for (int j = 0; j < instance.cardN; ++j)
								{
									if (Cplex.getValue(X[i][j]) >= 0.98) {
										ZZ1 = ZZ1 + Cplex.getValue(X[i][j]) * instance.C[i][j];

									}
								}

							float ZZ2;
							ZZ2 = 0;

							for (int i = 0; i < instance.cardN; ++i)
								for (int j = 0; j < instance.cardN; ++j)
								{
									if (Cplex.getValue(y[i][j]) >= 0.98) {
										ZZ2 = ZZ2 + Cplex.getValue(y[i][j]) * instance.C[i][j];
									}
								}

							//Grafico soluciones
							//Graficar_yed(argv[1], instance, Cplex, X, t, y, ZZ1, ZZ2);
						}

						/* Values to file */

						//Graficar(MIP_model, cplex_MIP, instance, X, t);
					}
					else
					{
						// error resolviendo la instancia, tomar acciones aca
					}
				}
				else
				{
					// error creando la instancia del modelo, tomar acciones aca
				}
				//cplex_MIP.exportModel("modelo2.lp"); // Exporta el modelo a un archivo.lp xxxx
				Model.end(); // eliminar el modelo
				Cplex.end(); // eliminar la instancia de CPLEX
				/* Tirar la cadena */
			}
			else
			{
				// error cargando datos, tomar medidas aca
			}
		}
		else
		{
			// error con el formato de la linea de comandos
			cout << "Error Entrada, use formato: [Ejecutable] [Instancia] [Rutas] [Soluciones]" << endl;
		}

	}
	catch (IloException &e)
	{
		env.out() << "ERROR: " << e << endl;
	}
	catch (...)
	{
		env.out() << "Excepcion no conocida" << endl;
	}
	env.end();
	return 0;
}

bool Load_data(IloEnv &env, char data_filename[255], Data &instance) // carga datos
{
	ifstream f_input(data_filename, ios::in); // abrir archivo de datos
	if (f_input)
	{
		f_input >> instance.cardN; // guardar cardN...

		// Symetric Instance

		instance.id = IloNumArray(env, instance.cardN);
		instance.x_cord = IloNumArray(env, instance.cardN);
		instance.y_cord = IloNumArray(env, instance.cardN);
		instance.SS = IloNumArray(env, instance.cardN);

		//Leo parametros
		for (int i = 0; i < instance.cardN; i++)
		{
			f_input >> instance.id[i] >> instance.x_cord[i] >> instance.y_cord[i];
		}

		instance.C = IloArray<IloNumArray>(env, instance.cardN);
		for (int i = 0; i < instance.cardN; i++)
		{
			instance.C[i] = IloNumArray(env, instance.cardN);
		}

		for (int i = 0; i < instance.cardN; i++)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				instance.C[i][j] = sqrt(pow(instance.x_cord[j] - instance.x_cord[i], 2) + pow(instance.y_cord[j] - instance.y_cord[i], 2));
			}
		}

		f_input >> instance.clusterN; // lee numero de cluster..

		vector<int> S1;
		vector<int> S2;
		vector<vector <int>> S3;
		IloNum S;
		int h = 0;
		int j = 0;

		for (int i = 0; i < instance.clusterN; i++)
		{
			//ingreso nodo 1(depot) al cluster 1
			if (i == 0)
			{
				S1.push_back(0);
				S3.push_back(S1);
			}
			S1.clear();
			S2.clear();
			S = 0;
			h = 0;
			while (S != -1)
			{
				f_input >> S;
				S1.push_back(S - 1);

				if (h != 0 && S == 1)
					S1.pop_back();
				h++;
			}

			if (S1.size() == 2)
				j++;

			if (S1.size() > 2)
			{
				for (int j = 0; j < S1.size() - 2; j++)
				{
					S2.push_back(S1[j + 1]);
				}
				S3.push_back(S2);
			}
		}

		//se le añade uno dado el nodo depósito
		if (j == 0)
			instance.clusterN++;

		instance.Cluster = IloArray<IloNumArray>(env, instance.clusterN);
		for (int i = 0; i < instance.Cluster.getSize(); i++)
		{
			instance.Cluster[i] = IloNumArray(env, S3[i].size());
			for (int j = 0; j < S3[i].size(); j++)
			{
				instance.Cluster[i][j] = S3[i][j];
			}
		}

		// En SS guardo el cluster al cual pertenece el nodo casillero
		for (int i = 0; i < instance.Cluster.getSize(); i++)
		{
			for (int j = 0; j < instance.Cluster[i].getSize(); j++)
			{
				instance.SS[instance.Cluster[i][j]] = i;
			}
		}

		// Mostrar datos
		/*
		cout << endl << "Cardinalidad: " << instance.cardN << endl;

		cout << endl << "id cord(X) cord(y):";

		for (int i = 0; i < instance.cardN; i++)
		{
			cout << endl << " " << instance.id[i] << " " << instance.x_cord[i] << " " << instance.y_cord[i];
		}

		cout << endl << "Matriz de costos:" << endl;
		for (int i = 0; i < instance.cardN; i++)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				cout << " " << instance.C[i][j];
			}
			cout << endl;
		}

		cout << endl << "Clusters: " << instance.clusterN << endl;

		cout << endl << "Matriz de clusters: " << endl;
		for (int i = 0; i < instance.Cluster.getSize(); i++)
		{
			for (int j = 0; j < instance.Cluster[i].getSize(); j++)
			{
				cout << " " << instance.Cluster[i][j];
			}
			cout << endl;
		}

		*/
		// Asymmetric Instance


		return true; // retornar el valor de que lo logro
	}
	else // en otro caso
	{
		cerr << "Problemas abriendo archivo " << data_filename << endl; // alegar que no se pudo
		return false; // retornar la respuesta que no se pudo
	}
}

bool Create_MIP(Data &instance, IloModel &Model, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y, IloIntVarArray &t) // crea la instancia del modelo
{
	IloEnv env = Model.getEnv(); // cachar el entorno actual
	char varName[50]; // arreglo de caracteres para ponerle nombre a las variables    OJO depende del tamaño de la red

	//X_ij
	X = IloArray<IloBoolVarArray>(env, instance.cardN); // X es una matriz 2D de variables enteras de cardN filas

	for (int i = 0; i < instance.cardN; ++i) // cada fila
	{
		X[i] = IloBoolVarArray(env, instance.cardN); // es un arreglo de variables binarias de dimension cardN
		for (int j = 0; j < instance.cardN; ++j)
		{
			sprintf(varName, "X[%d,%d]", i, j); // Ponerle nombre a cada variable
			X[i][j].setName(varName);
			if (i == j)
				X[i][j].setBounds(0, 0);
		}
		Model.add(X[i]); // Agregarla al modelo
	}

	//y_ij
	y = IloArray<IloBoolVarArray>(env, instance.cardN);

	for (int i = 0; i < instance.cardN; ++i) // cada fila
	{
		y[i] = IloBoolVarArray(env, instance.cardN); // es un arreglo de variables binarias de dimension cardN
		for (int j = 0; j < instance.cardN; ++j)
		{
			sprintf(varName, "y[%d,%d]", i, j); // Ponerle nombre a cada variable
			y[i][j].setName(varName);
		}
		Model.add(y[i]); // Agregarla al modelo
	}

	//fijo en valor 0 las asignaciones a nodos que no son del cluster
	for (int p = 0; p < instance.clusterN; p++)
	{
		for (int q = 0; q < instance.clusterN; q++)
		{
			if (q != p)
			{
				for (int i = 0; i < instance.Cluster[p].getSize(); i++)
				{
					for (int j = 0; j < instance.Cluster[q].getSize(); j++)
					{
						y[instance.Cluster[p][i]][instance.Cluster[q][j]].setBounds(0, 0);
					}
				}
			}
		}
	}

	//t_i
	t = IloIntVarArray(env, instance.clusterN, 1, instance.clusterN - 1); //Definir como double ????

	for (int i = 0; i < instance.clusterN; i++)
	{
		sprintf(varName, "t[%d]", i);
	}
	Model.add(t);

	cout << endl << "Definicion de variables ok" << endl;

	//Definir FO del problema
	IloExpr Obj1(env);                                   // expresion para el objetivo
	IloExpr Obj2(env);
	for (int i = 0; i < instance.cardN; ++i)
		for (int j = 0; j < instance.cardN; ++j)
			if (i != j) {
				Obj1 += (instance.C[i][j] * X[i][j]);
				Obj2 += (instance.C[i][j] * y[i][j]);
			}
	Model.add(IloMinimize(env, Obj1 + Obj2)); // Z = min FO
	Obj1.end();                            // cerrar la expresion
	Obj2.end();                            // cerrar la expresion

	/**********************Restricciones del Modelo de Asignación (Por cluster)**********************/

	// Restriccion (2) sale de cluster Np
	int Cont = 0;
	for (int p = 0; p < instance.clusterN; ++p)
	{
		IloExpr Out_deeg(env);
		for (int i = 0; i < instance.Cluster[p].getSize(); ++i)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				Cont = 0;
				for (int h = 0; h < instance.Cluster[p].getSize(); ++h)   //veo si el nodo j esta en el cluster p
				{
					if (j == instance.Cluster[p][h])
					{
						Cont = 1;
						break;
					}
				}
				if (Cont == 0)
					Out_deeg += X[instance.Cluster[p][i]][j];
			}
		}
		Model.add(Out_deeg == 1); // sum{j} X_ij = 1, para todo i
		Out_deeg.end();
	}

	// Restriccion (3) entra a p
	for (int p = 0; p < instance.clusterN; ++p)
	{
		IloExpr In_deeg(env);
		for (int i = 0; i < instance.Cluster[p].getSize(); ++i)
		{
			for (int j = 0; j < instance.cardN; j++)
			{
				Cont = 0;
				for (int h = 0; h < instance.Cluster[p].getSize(); ++h)   //veo si el nodo j esta en el cluster p
				{
					if (j == instance.Cluster[p][h])
					{
						Cont = 1;
						break;
					}
				}
				if (Cont == 0)
					In_deeg += X[j][instance.Cluster[p][i]];
			}
		}
		Model.add(In_deeg == 1); // sum{j} X_ij = 1, para todo i
		In_deeg.end();
	}

	// Restriccion (4) balance en todos los nodos del tour
	for (int j = 0; j < instance.cardN; ++j)
	{
		IloExpr Out_deeg(env);
		IloExpr In_deeg(env);
		for (int i = 0; i < instance.cardN; ++i)
		{
			if (i != j) {
				Out_deeg += X[i][j];
				In_deeg += X[j][i];
			}
		}
		Model.add(Out_deeg == In_deeg); // sum{i} X_ij = sum{h} X_jh para todo j
		Out_deeg.end();
		In_deeg.end();
	}

	// Restriccion (5) el tour tiene tantos arcos como cluster hay. Co esta restriccion se impone que el tour visita un solo nodo de cada cluster 
	/*
	IloExpr arcs(env);
	for (int i = 0; i < instance.cardN; ++i)
	{
		for (int j = 0; j < instance.cardN; ++j)
		{
			if (i != j)
				arcs += X[i][j];
		}
	}
	Model.add(arcs >= instance.clusterN); //
	arcs.end();
	*/
	// Restriccion (6) todos los nodos de un cluster son asignados al paradero dentro del cluster
	for (int p = 0; p < instance.clusterN; ++p)
	{
		for (int i = 0; i < instance.Cluster[p].getSize(); ++i)
		{
			IloExpr Assigns(env);
			for (int j = 0; j < instance.Cluster[p].getSize(); ++j)
			{
				Assigns += y[instance.Cluster[p][i]][instance.Cluster[p][j]];
			}
			Model.add(Assigns == 1); // 
			Assigns.end();
		}
	}

	// Restriccion (7) captura el paradero dentro del cluster
	for (int p = 0; p < instance.clusterN; ++p)
	{
		for (int j = 0; j < instance.Cluster[p].getSize(); ++j)
		{
			IloExpr arcs(env);
			for (int i = 0; i < instance.cardN; ++i)
			{
				if (i != instance.Cluster[p][j])
				{
					arcs += X[i][instance.Cluster[p][j]];
				}
			}
			Model.add(arcs == y[instance.Cluster[p][j]][instance.Cluster[p][j]]); // 
			arcs.end();
		}
	}

	// Restriccion (8)	   
	for (int p = 0; p < instance.clusterN; ++p)
	{
		for (int i = 0; i < instance.Cluster[p].getSize(); ++i)
		{
			for (int j = 0; j < instance.Cluster[p].getSize(); ++j)
			{
				if (i != j)
				{
					Model.add(y[instance.Cluster[p][i]][instance.Cluster[p][j]] <= y[instance.Cluster[p][j]][instance.Cluster[p][j]]); // 
				}
			}
		}
	}

	cout << endl << "Restricciones de asignacion OK " << endl;

	/**************Restriciones de Miller-Tucker-Zemlin (por cluster)***************/

	//Res (12): Primer nodo
	t[0].setBounds(0, 0);

	//Res (13): definicion de variable
	for (int p = 1; p < instance.clusterN; p++)
	{
		Model.add(1 <= t[p] <= instance.clusterN - 1);
	}

	//Res (14): Eliminacion de subtour
	for (int p = 0; p < instance.clusterN; p++)
	{
		for (int q = 1; q < instance.clusterN; q++)
		{
			if (q != p)
			{
				IloExpr In_deeg(env);
				for (int i = 0; i < instance.Cluster[p].getSize(); i++)
				{
					for (int j = 0; j < instance.Cluster[q].getSize(); j++)
					{
						In_deeg += X[instance.Cluster[p][i]][instance.Cluster[q][j]];
					}
				}
				//cout << endl << " Salto " << Arcs_cluster << endl;
				Model.add(t[q] >= t[p] + In_deeg + (instance.clusterN - 1) * (In_deeg - 1));
				In_deeg.end();
			}
		}
	}

	cout << endl << "Restricciones de Miller-Tucker-Zemlin por cluster OK " << endl;

	return true;
}

bool MIP_solve(IloModel &Model, IloCplex &Cplex, Data &instance, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y, IloIntVarArray &t)
{

	//cout << "MIP_solve" << endl;
	IloEnv env = Model.getEnv(); // cachar el entorno


	//Cplex.setParam(cplex_MIP.PreInd, 0); // apagar presolve
	//Cplex.setOut(env.getNullStream()); // apaga reporte de CPLEX por pantalla
	//cplex_MIP.setWarning(env.getNullStream()); // apaga alegatos precautorios de CPLEX
	//cplex_MIP.setParam(IloCplex::ClockType, 1); // 1: cpu time // 2: wall clock time // 0: defaul (wall clock time)
	Cplex.setParam(IloCplex::ClockType, 1); // 1: cpu time // 2: wall clock time // 0: defaul (wall clock time)
	Cplex.setParam(IloCplex::TiLim, 3600); // limite de tiempo en segundos
	Cplex.setParam(IloCplex::ParallelMode, IloCplex::Deterministic);
	Cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
	/*
	cplex_MIP.setParam(IloCplex::MIPSearch, IloCplex::Traditional);//B&B tradicional (no busqueda dinamica)
	cplex_MIP.setParam(IloCplex::ParallelMode, IloCplex::Deterministic);//B&B deterministico(no oportunista)

	cplex_MIP.setParam(IloCplex::MIPDisplay, 2);
	cplex_MIP.setParam(IloCplex::Threads, 1);
	cplex_MIP.setParam(IloCplex::MIPSearch, 1);
	cplex_MIP.setParam(IloCplex::PreInd, 0);
	cplex_MIP.setParam(IloCplex::CutsFactor, 0);
	cplex_MIP.setParam(IloCplex::Cliques, -1);
	cplex_MIP.setParam(IloCplex::Covers, -1);
	cplex_MIP.setParam(IloCplex::DisjCuts, -1);
	cplex_MIP.setParam(IloCplex::FlowCovers, -1);
	cplex_MIP.setParam(IloCplex::FlowPaths, -1);
	cplex_MIP.setParam(IloCplex::FracCuts, -1);
	cplex_MIP.setParam(IloCplex::GUBCovers, -1);
	cplex_MIP.setParam(IloCplex::ImplBd, -1);
	cplex_MIP.setParam(IloCplex::MIRCuts, -1);
	cplex_MIP.setParam(IloCplex::ZeroHalfCuts, -1);
	cplex_MIP.setParam(IloCplex::EachCutLim, 0);
	cplex_MIP.setParam(IloCplex::CutPass, -1);
	cplex_MIP.setParam(IloCplex::HeurFreq, -1);
	*/

	//aplico cortes enteros
	eps_rhs = Cplex.getParam(IloCplex::EpRHS);

	FacilityCallback fcCallback(X, y, instance);

	CPXLONG contextMask = 0;
	contextMask |= IloCplex::Callback::Context::Id::Candidate;

	// If contextMask is not zero we add the callback.
	if (contextMask != 0)
		Cplex.use(&fcCallback, contextMask);

	if (Cplex.solve()) // si se pudo resolver avisar
		return true;

	else // en otro caso, tambien avisar
		return false;
}


void Guardar(char *nombre, char resultados[50], IloCplex &Cplex, double &time, IloArray<IloBoolVarArray> &X, IloArray<IloBoolVarArray> &y)
{
	string instancia; //guarda nombre de la instancia
	int i = 0;
	while (nombre[i] != '.')
	{
		instancia.push_back(nombre[i]);
		i++;
	}

	time = Cplex.getTime();

	fstream salida2(resultados, ios::app);
	if (salida2.is_open())
	{
		salida2 << instancia << "\t" << Cplex.getObjValue() << "\t" << Cplex.getMIPRelativeGap() << "\t" << Cplex.getNnodes() << "\t" << time << endl;
	}
	else {
		cerr << "Problemas guardando soluciones " << endl;
	}

	salida2.close();
}

void Graficar_yed(char *nombre, Data &instance, IloCplex &Cplex, IloArray<IloBoolVarArray> &X, IloIntVarArray &t, IloArray<IloBoolVarArray> &y, float &ZZ1, float &ZZ2)
{
	IloEnv env = Cplex.getEnv();
	//paso valores de variables

	IloArray<IloNumArray> X_Sol(env, instance.cardN);
	for (int i = 0; i < instance.cardN; ++i)
	{
		X_Sol[i] = IloNumArray(env, instance.cardN);
		Cplex.getValues(X_Sol[i], X[i]);
	}

	IloArray<IloNumArray> Y_Sol(env, instance.cardN);
	for (int i = 0; i < instance.cardN; ++i)
	{
		Y_Sol[i] = IloNumArray(env, instance.cardN);
		Cplex.getValues(Y_Sol[i], y[i]);
	}

	//para crear gml con nomobre de la instancia
	string instancia; //guarda nombre de la instancia
	int i = 0;
	while (nombre[i] != '.')
	{
		instancia.push_back(nombre[i]);
		i++;
	}

	string archivo_salida = instancia + ".gml";

	ofstream gml(archivo_salida, ios::out);
	double Factorx = 2.2;
	double Factory = 2.6;
	double Factor2x = 3.7;
	double Factor2y = 3.5;
	int Depot = 1;
	gml << " graph [hierachic 1 directed 1 " << endl;


	//calculo centroide y lo amplio para incluir nodos
	vector <int> Cx, Cy, H, W;
	int px, py;
	int ax, bx;
	int ay, by;
	int width, high = 0.0;

	for (int i = 0; i < instance.Cluster.getSize(); i++)
	{
		ax = 0;
		bx = 0;
		ay = 0;
		by = 0;
		px = 0;
		py = 0;
		width = 0;
		high = 0;
		//cout << endl << "cluster: " << i << " Nodos: " << endl;
		for (int j = 0; j < instance.Cluster[i].getSize(); j++)
		{
			//cout << " " << instance.Cluster[i][j];
			if (j == 0)
			{
				ax = instance.x_cord[instance.Cluster[i][j]];
				bx = instance.x_cord[instance.Cluster[i][j]];
				ay = instance.y_cord[instance.Cluster[i][j]];
				by = instance.y_cord[instance.Cluster[i][j]];
			}

			//cout << endl << " Coordenada (x): " << instance.xcoord[instance.Cluster[i][j]] << " Coordenada (y): " << instance.ycoord[instance.Cluster[i][j]];
			if (j > 0 && instance.x_cord[instance.Cluster[i][j]] > ax)
				ax = instance.x_cord[instance.Cluster[i][j]];

			if (j > 0 && instance.x_cord[instance.Cluster[i][j]] < bx)
				bx = instance.x_cord[instance.Cluster[i][j]];

			if (j > 0 && instance.y_cord[instance.Cluster[i][j]] > ay)
				ay = instance.y_cord[instance.Cluster[i][j]];

			if (j > 0 && instance.y_cord[instance.Cluster[i][j]] < by)
				by = instance.y_cord[instance.Cluster[i][j]];
		}
		//calculo centroides
		if (instance.Cluster[i].getSize() > 1)
		{
			px = ((ax - bx) / 2) + bx;
			py = ((ay - by) / 2) + by;
		}

		if (instance.Cluster[i].getSize() == 1)
		{
			px = ax;
			py = ay;
		}

		//agrego centroides
		Cx.push_back(px);
		Cy.push_back(py);

		//calculo ancho y alto
		px = 0;
		py = 0;

		if (instance.Cluster[i].getSize() > 1)
		{
			width = (ax - bx + 30) * Factor2x;
			high = (ay - by + 20) * Factor2y;
		}

		if (instance.Cluster[i].getSize() == 1)
		{
			width = 45 * Factorx;
			high = 45 * Factory;
		}

		//agrego alto y ancho a vectores
		W.push_back(width);
		H.push_back(high);
	}

	//Grafico clusters
	for (int i = 0; i < instance.clusterN; i++)
	{
		gml << "node [ id " << instance.cardN + i + 2 << " graphics [ x " << Cx[i] * Factorx << " y " << Cy[i] * (-1)* Factory
			<< " w " << W[i] << " h " << H[i] << " type \"ellipse\" fill \"#ffff99\"] LabelGraphics [text \"\" fontSize 30 ] ]" << endl;
	}

	//grafico nodo deposito
	for (int i = 0; i < Depot; ++i)
	{
		gml << "node [ id " << i << " graphics [ x " << instance.x_cord[i] * Factorx
			<< " y " << instance.y_cord[i] * (-1) * Factory
			<< " w 60 h 60 type \"ellipse\" fill \"#ff9900\"] LabelGraphics [text \""
			<< 'D' << "\" fontSize 30 ] ]" << endl;
	}

	//grafico el resto de los nodos
	for (int i = 1; i < instance.cardN; ++i)
	{
		gml << "node [ id " << i << " graphics [ x " << instance.x_cord[i] * Factorx
			<< " y " << instance.y_cord[i] * (-1) * Factory
			<< " w 60 h 60 type \"ellipse\" fill \"#ccffff\"] LabelGraphics [text \""
			<< i + 1 << "\" fontSize 30 ] ]" << endl;
	}


	//Grafico los Arcos entre Clusters
	for (int i = 0; i < instance.cardN; i++)
	{
		for (int j = 0; j < instance.cardN; j++)
		{
			if (j != i && X_Sol[i][j] > 1e-6 && X_Sol[i][j] <= 1.0 + 1e-6)
			{
				gml << "edge [source " << i << " target " << j
					<< " graphics [ width 3 style \"solid\"	fill \"#000000\" targetArrow \"standard\" ] LabelGraphics [ text \"\" fontSize 25 ] ]" << endl;
			}
		}
	}

	//Gafico los Arcos dentro del Cluster
	for (int i = 0; i < instance.cardN; i++)
	{
		for (int j = 0; j < instance.cardN; j++)
		{
			if (j != i && Y_Sol[i][j] > 1e-6 && Y_Sol[i][j] <= 1.0 + 1e-6) //grafico los fraccionarios
			{

				gml << "edge [source " << i << " target " << j
					<< " graphics [ width 4 style \"dashed\"	fill \"#ff6600\" targetArrow \"standard\" ] LabelGraphics [ text \"\" fontSize 25 ] ]" << endl;
			}
		}
	}

	//grafico las Funciones objetivos

	int X_cord = 0;
	int Y_cord = 0;
	for (int i = 0; i < instance.cardN; i++)
	{
		if (instance.x_cord[i] < X_cord)
			X_cord = instance.x_cord[i];

		if (instance.y_cord[i] < Y_cord)
			Y_cord = instance.y_cord[i];
	}

	gml << " node [ id 79 graphics [ x  " << (X_cord * Factorx) + 200 << "  y " << (Y_cord * (-1) * Factory) << " w 350 h 60 type \"rectangle\" fill \"#ccffff\"] LabelGraphics [text \" Z1 = " << ZZ1 << "\" fontSize 30 ] ]" << endl;

	gml << " node [ id 79 graphics [ x  " << (X_cord * Factorx) + 200 << "  y " << (Y_cord * (-1) * Factory) - 90 << " w 350 h 60 type \"rectangle\" fill \"#ccffff\"] LabelGraphics [text \" Z2 = " << ZZ2 << "\" fontSize 30 ] ]" << endl;

	gml << "]" << endl;
	gml.close();
}