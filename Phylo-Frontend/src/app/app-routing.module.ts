import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';

// General imports
import { HomeComponent } from './layout/home/home.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';

const routes: Routes = [
  {path: 'home', component: HomeComponent},
  {path: 'loadGof', component: GofHolderComponent},
  {
      path: '**',
    redirectTo: '/home',
    pathMatch: 'full'
  }
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
