/*
 * CalculsRigidBody.cpp :
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
#include "ObjetSimuleRigidBody.h"
#include "Viewer.h"

using namespace std;


/*
 * Calcul de la masse de l objet rigide
 */
void ObjetSimuleRigidBody::CalculMasse()
{
	_Mass = 0.0;
	for (int i = 0; i < P.size(); i++)
	{
		_Mass += M[i];
		_BaryCentre = _BaryCentre + P[i] * M[i];
	}
	_BaryCentre = _BaryCentre / _Mass;
}

/*
 * Calcul du tenseur d inertie de l objet rigide - - partie constante Ibody
 * et remplissage tableau roi (position particules dans repere objet)
 */
void ObjetSimuleRigidBody::CalculIBody()
{
	_Ibody = Matrix::NullMatrix();
	for (int i = 0; i < P.size(); i++)
	{
		_ROi.push_back(P[i]);

		float r0iTr0i = (_ROi[i].x * _ROi[i].x + _ROi[i].y * _ROi[i].y + _ROi[i].z * _ROi[i].z);
		Matrix r0ir0iT = MultiplyTransposedAndOriginal(_ROi[i]);
		_Ibody = (Matrix::UnitMatrix() * r0iTr0i - r0ir0iT) * M[i];
	}

	_IbodyInv = _Ibody.InverseConst();
}

/*
 * Calcul de l etat de l objet rigide.
 */
void ObjetSimuleRigidBody::CalculStateX()
{
	_InertieTenseurInv = _Rotation * _IbodyInv * _Rotation.TransposeConst();
	_VitesseAngulaire = _InertieTenseurInv * _MomentCinetique;
}

/*
 * Calcul de la derivee de l etat : d/dt X(t).
 */
void ObjetSimuleRigidBody::CalculDeriveeStateX(Vector gravite)
{
	_Vitesse = _QuantiteMouvement / _Mass;
	_RotationDerivee = (_VitesseAngulaire.x == 0.0 && _VitesseAngulaire.y == 0.0 && _VitesseAngulaire.z == 0.0) ? Matrix::NullMatrix() : StarMatrix(_VitesseAngulaire) * _Rotation;
	_Force = gravite * _Mass;

	for (int i = 0; i < _Nb_Sommets; i++)
		_Torque = _Torque + cross(P[i] - _BaryCentre, F[i]);
}

/**
 * Schema integration pour obbtenir X(t+dt) en fct de X(t) et d/dt X(t)
 */
void ObjetSimuleRigidBody::Solve(float visco)
{
	_Position = _Position + _Vitesse * _delta_t;
	_Rotation = _Rotation + _RotationDerivee * _delta_t;
	_QuantiteMouvement = _QuantiteMouvement + _Force * _delta_t;
	_MomentCinetique = _MomentCinetique + _Torque * _delta_t;

	for (int i = 0; i < _Nb_Sommets; i++)
		P[i] = _Rotation * _ROi[i] + _Position * visco;

}//void


void ObjetSimuleRigidBody::CollisionPlan()
{

}

/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleRigidBody::CollisionPlan(float x, float y, float z)
{
	for (int i = 0; i < V.size(); i++)
	{
		if (P[i].y <= y && _QuantiteMouvement.y < 0.0)
		{
			_QuantiteMouvement = -_QuantiteMouvement * 0.5f;

			if (_MomentCinetique.x < 0.01 || _MomentCinetique.y < 0.01 || _MomentCinetique.z < 0.01)
				_MomentCinetique = Vector(0.0, 0.0, 0.0);

			_MomentCinetique = -_MomentCinetique * 0.5f;
		}
	}
}// void