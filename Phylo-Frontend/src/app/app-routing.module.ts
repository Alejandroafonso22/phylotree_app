import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

// General imports
import { HomeComponent } from './layout/home/home.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { SpeciesdbComponent } from './speciesdb/speciesdb.component';
import {Molecule3DComponent} from './molecule3-d/molecule3-d.component';
import {NcbiinsightsComponent} from './ncbiinsights/ncbiinsights.component';

const routes: Routes = [
  {path: 'about', component: SpeciesdbComponent},
  {path: 'molecule3D', component: Molecule3DComponent},
  {path: 'insights', component: NcbiinsightsComponent},
  {path: 'home', component: HomeComponent},
  {path: 'loadGof', component: GofHolderComponent},
  {path: '**', component: NotFoundComponent},
  {path: '', redirectTo: 'home', pathMatch: 'full'},


];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
