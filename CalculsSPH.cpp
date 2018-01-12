/*
 * CalculsSPH.cpp :
 * Copyright (C) 2017 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 * Utilisation du code :
 * https://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
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


 /** \file CalculsSPH.cpp
 Programme calculant pour  un fluide son etat au pas de temps suivant (methode d 'Euler semi-implicite) en utilisant la methode SPH (Smoothed Particles Hydrodynamics):
 principales fonctions de calculs.
 \brief Fonctions de calculs pour un fluide avec methode SPH.
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleSPH.h"
#include "Viewer.h"

using namespace std;

/**
 * Calcul des densites des particules.
 * Formule :
 *  \rho_i = \frac{4m}{\pi h^8} \sum_{j \in N_i} (h^2 - r^2)^3.
 */
void ObjetSimuleSPH::CalculDensite()
{
	for (int i = 0; i < rho.size(); i++)
	{
		rho[i] = 0;
	}

	for (int i = 0; i < rho.size(); i++)
	{
		float _rho = (4 * M[i]) / (M_PI * pow(h, 8));

		float sum = 0;
		for (int j = i + 1; j < rho.size(); j++)
		{
			float r = distance(Point(P[i]), Point(P[j]));

			float d = (h * h) - (r *r);
			if (d > 0)
				sum += d * d * d;

			rho[i] += _rho * sum;
			rho[j] += _rho * sum;
		}
	}

}//void

/**
 * Calcul des forces d interaction entre particules.
  * Attention - Calcul direct de fij / rho_i
 */
void ObjetSimuleSPH::CalculInteraction(float visco)
{
	for (int i = 0; i < rho.size(); i++)
	{
		for (int j = i + 1; j < rho.size(); j++)
		{
			float r = distance(Point(P[i]), Point(P[j]));
			float d = (h * h) - (r *r);
			if (d > 0)
			{
				float qij = r / h;
				Vector rij = P[i] - P[j];
				Vector vij = V[i] - V[j];

				F[i] = (M[j] / (M_PI * pow(h, 4) * rho[j])) * (1 - qij) *
					((15 * bulk) * (rho[i] + rho[j] - 2 * rho0) * ((1 - qij) / qij) * rij - (40 * visco * vij));

				F[i] = F[i] / rho[i];
				F[j] = -F[i];
			}
		}
	}
}//void

/**
 * Gestion des collisions.
 * Notre condition aux limites correspond à une frontière inélastique
 * avec un coefficient de restitution spécifié < 1. Quand
 * une particule frappe une barrière verticale ([[which = 0]]) ou horizontale
 * barrier ([[which = 1]]), nous le traitons avec [[damp_reflect]].
 * Cela réduit la distance totale parcourue en fonction du temps écoulé depuis
 * la collision reflète, amortit les vitesses, et reflète
 * quels que soient les composants de la solution à refléter.
 */
void ObjetSimuleSPH::damp_reflect(int frontiere, float barrier, int indice_part)
{
	/// frontiere : indique quelle frontiere (x, y, z) du domaine est concernee
	const float DAMP = 0.75;


	float tbounce = (P[indice_part](frontiere) - barrier) / V[indice_part](frontiere);

	P[indice_part].x = P[indice_part].x - (V[indice_part].x * (1 - DAMP) * tbounce);
	P[indice_part].y = P[indice_part].y - (V[indice_part].y * (1 - DAMP) * tbounce);
	P[indice_part].z = P[indice_part].z - (V[indice_part].z * (1 - DAMP) * tbounce);


	P[indice_part](frontiere) = 2 * barrier - P[indice_part](frontiere);
	V[indice_part](frontiere) = -V[indice_part](frontiere);	V[indice_part](0) = V[indice_part](0) * DAMP;	V[indice_part](1) = V[indice_part](1) * DAMP;	V[indice_part](2) = V[indice_part](2) * DAMP;
}

/**
 * Gestion des collisions.
 * Pour chacune des particules nous verifions la reflection
 * avec chacun des 4 murs du domaine.
 */
void ObjetSimuleSPH::CollisionPlan(float x, float y, float z)
{
}