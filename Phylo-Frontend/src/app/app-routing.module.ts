import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

// General imports
import { HomeComponent } from './layout/home/home.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { SpeciesComponent } from './species/species.component';
import { SpeciesdbComponent } from './speciesdb/speciesdb.component';

const routes: Routes = [
  {path: 'speciesdb', component: SpeciesdbComponent},
  {path: 'home', component: HomeComponent},
  {path: 'species', component: SpeciesComponent},
  {path: 'loadGof', component: GofHolderComponent},
  {path: '**', component: NotFoundComponent},
  {path: '', redirectTo: 'home', pathMatch: 'full'},


];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
