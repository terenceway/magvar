/*
 * https://www.wayforward.net/magvar/av_mag.h
 * Copyright(c) 2005, 2008, 2022 Wayforward LLC
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
#ifndef AV_MAG_H
#define AV_MAG_H

#ifndef PI
# define PI (3.14159265358979323f)
#endif

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * Calculate the magnetic variation given a Julian date,
	 * the latitude and longitude, and the elevation above MSL.
	 *
         * @returns magnetic variation, in radians.
	 *          Negative values are EAST deviation "east is least,"
	 *          positive values are WEST deviation.  This is
	 *          different from the NGDC software!
         *
	 * @param sdate  floating point julian date (for example,
	 *               2004.92 is in December, 2004)
	 * @param latitude latitude in radians (-pi/2..pi/2)
	 * @param longitude longitude in radians (-pi..pi)
	 * @param elev elevation in meters
	 */
	float
	magnetic_variation(float sdate, float latitude,
			   float longitude, float elev);

#ifdef __cplusplus
}
#endif

#endif /* AV_MAG_H */

/*
 * Overrides for Emacs so that we follow Linux's tabbing style.
 * Emacs will notice this stuff at the end of the file and automatically
 * adjust the settings for this buffer only.  This must remain at the
 * end of the file.
 * ----------------------------------------------------------------------
 * Local Variables:
 * c-file-style: "linux"
 * End:
 */
