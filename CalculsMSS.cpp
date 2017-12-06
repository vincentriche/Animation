/*
 * CalculsMSS.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

 /** \file CalculsMSS.cpp
 Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant
  (methode d 'Euler semi-implicite) : principales fonctions de calculs.
 \brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;

/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext

	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext

	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration

	std::vector<Ressort*> ressorts_List = _SytemeMasseRessort->GetRessortList();
	for (int i = 0; i < ressorts_List.size(); i++)
	{
		Ressort* ressort = ressorts_List[i];
		Particule* particleA = ressort->GetParticuleA();
		Particule* particleB = ressort->GetParticuleB();

		float distanceAB = distance(Point(particleA->GetPosition()), Point(particleB->GetPosition()));
		if (distanceAB > ressort->GetLrepos() * 2.0f)
		{
			_SytemeMasseRessort->GetRessortList().erase(_SytemeMasseRessort->GetRessortList().begin() + i);
			continue;
		}

		float cnst_raideur = ressort->GetRaideur();
		float cnst_amortissement = ressort->GetAmortissement();

		Vector vectorUij = normalize(P[particleB->GetId()] - P[particleA->GetId()]);

		Vector forcesE = cnst_raideur * (distance(Point(P[particleB->GetId()]), Point(P[particleA->GetId()]))
			- ressort->GetLrepos()) * vectorUij;

		Vector forcesV = (cnst_amortissement * dot(V[particleB->GetId()] - V[particleA->GetId()], vectorUij))
			* vectorUij;

		Vector wind = Vector(0.0, 0.0, 0.0);

		Vector somme_forces = forcesE + forcesV + wind;

		F[particleA->GetId()] = F[particleA->GetId()] + somme_forces;
		F[particleB->GetId()] = F[particleB->GetId()] - somme_forces;
	}

}//void

 /**
 * Gestion des collisions avec le sol - plan (x,y,z).
 */
void ObjetSimuleMSS::CollisionPlan(float x, float y, float z)
{

	for (int i = 0; i < V.size(); i++)
	{
		if (P[i].y < y)
		{
			P[i] = Vector(P[i].x, y, P[i].z);
			V[i] = Vector(0, 0, 0);
		}
	}
}// void

 /**
 * Gestion des collisions avec une sphère.
 */
void ObjetSimuleMSS::CollisionSphere(Point p, float r)
{
	for (int i = 0; i < P.size(); i++)
	{
		Vector direction = normalize(Vector(P[i]) - Vector(p));
		if (distance(Point(P[i]), p) < r)
		{
			P[i] = Vector(p) + Vector(direction) * r;
			V[i] = Vector(0, 0, 0);
		}
	}
}// void

