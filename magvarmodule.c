/*
 * https://www.wayforward.net/magvar/magvarmodule.c
 * Copyright(c) 2005, 2008, 2022 Wayforward LLC
 *
 * Python3 C module.
 *
 * ---
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "av_mag.h"

static PyObject *
magvar_magnetic_variation(PyObject *self, PyObject *args) {
	static const char ARGUMENTS[] = "ffff";

	float date;
	float latitude;
	float longitude;
	float elevation;

	float magvar;

	if (!PyArg_ParseTuple(args, ARGUMENTS, &date, &latitude,
			      &longitude, &elevation))
		return NULL;

	magvar = magnetic_variation(date, latitude, longitude, elevation);

	return PyFloat_FromDouble(magvar);
}

static PyMethodDef magvar_methods[] = {
	{"magnetic_variation", magvar_magnetic_variation, METH_VARARGS,
	  "Return magnetic variation, in radians, given date (julian),"
	  " latitude (radians), longitude (radians), and elevation (meters)"},

	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef magvarmodule = {
	PyModuleDef_HEAD_INIT,
	"magvar",
	"Geomagnetic variation (declination), from the WMM COF data",
	-1,
	magvar_methods
};

PyMODINIT_FUNC
PyInit_magvar(void) {
	return PyModule_Create(&magvarmodule);
}
