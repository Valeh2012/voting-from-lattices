/**
 *
 * This file is part of FV-NFLlib
 *
 * Copyright (C) 2016  CryptoExperts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */

#include <cstddef>
#include <vector>

#include <gmpxx.h>
#include <nfl.hpp>

/// include the FV homomorphic encryption library



int main() {
	using P = nfl::poly_from_modulus<uint64_t, 2048, nfl::params<uint64_t>::kModulusBitsize>;
	P a = nfl::uniform();
	P b = nfl::uniform();
	P c = a + b;
	std::cout << c << std::endl;
	return 0;
}
