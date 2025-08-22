# Copyright © 2024, Sergio Callegari
# All rights reserved.

# This file is part of PyDSM.

# PyDSM is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyDSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

def pytest_configure(config):
    config.addinivalue_line("python_files", "benchmark_*.py")
    config.addinivalue_line("python_functions", "benchmark_*")
    config.addinivalue_line("python_classes", "Benchmark_*")
    config.addinivalue_line("markers",
        "slow: Tests that are very slow.")
