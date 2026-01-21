%define _topdir       /tmp/cbwavesrpm
%define name          cbwaves
%define version       1.2.0
%define release       2
%define buildroot     %{_topdir}/%{name}-%{version}-%{release}-root

BuildRoot: %{buildroot}
Summary:   %{name} - a gravitational wave template generator
Name:      %{name}
Version:   %{version}
Release:   %{release}
License:   GPL
Group:     Applications/Gravity
Source:    %{name}-%{version}-%{release}.tgz
URL:       http://virgo.rmki.kfki.hu/
Vendor:    Virgo Group of Wigner Research Centre for Physics
Packager:  Gergely Debreczeni <Debreczeni.Gergely@wigner.mta.hu>
Requires:  boost-devel

%description
The cbwaves software calculates the gravitational waveforms
emitted by general configuration (not-aligned spinning,
eccentric) BNSs or BBHs. It obtains the time domain waveforms 
by directly integrating the equation of motion of the bodies.

%prep
tar -xzvf ../SOURCES/%{name}-%{version}-%{release}.tgz 
cd %{name}-%{version}-%{release}

%build
cd %{name}-%{version}-%{release}/src
make

%install
mkdir -p %{_topdir}/%{name}-%{version}-%{release}-root
cd %{name}-%{version}-%{release}/src
make install prefix=%{buildroot}/usr/local/cbwaves

%files
/usr/local/cbwaves/

